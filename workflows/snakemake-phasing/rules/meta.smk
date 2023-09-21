import pandas
from snakemake.utils import validate
from snakemake.utils import min_version
import pickle
import typing
import os
from pathlib import Path
from snakemake.logging import logger
import pysam

###################
# Snakemake Setup #
###################
min_version("6.1.0")
report: "../report/phasing.rst"

configfile: "config.yml"
validate(config, schema="../schemas/config.schema.yml")

###############
# Input Files #
###############
if 'id-mapping' in config and 'inova' in config['id-mapping']:
    inova_id_match = pandas.read_csv(
        config['id-mapping']['inova'], skiprows=1, names=["usi", "rid"]
    ).set_index("usi", drop=False)
else:
    inova_id_match = pandas.DataFrame([], columns=["usi", "rid"])

if 'id-mapping' in config and 'radar' in config['id-mapping']:
    usi2rid = pandas.read_csv(
        config['id-mapping']['radar'], sep="\t", names=["usi", "rid"]
    ).set_index("usi", drop=False)
else:
    usi2rid = pandas.DataFrame([], columns=["usi", "rid"])

usi2rid = pandas.concat([usi2rid, inova_id_match], sort=True)

files = pandas.read_csv(
    config["files"],
    sep="\t",
    header=0,
).set_index("sample_id", drop=False)
validate(files, schema="../schemas/files.schema.yml")

def get_sample(sample: str, usi2rid=usi2rid):
    if sample in usi2rid.index:
        return usi2rid.loc[sample].rid
    return sample

def find_files(search_dir, suffix=".bam", index_suffix=".bai"):
    base = Path(search_dir)
    found_bams = []
    found_indexes = []
    for child in base.iterdir():
        if child.is_file():
            if child.suffix == suffix:
                found_bams.append(child)
            elif child.suffix == index_suffix:
                found_indexes.append(child)
        elif child.is_dir():
            other_bams, other_indexes = find_files(str(child),
                                                   suffix=suffix,
                                                   index_suffix=index_suffix)
            found_bams = [ *found_bams, *other_bams ]
            found_indexes = [ *found_indexes, *other_indexes ]
    return found_bams, found_indexes

sample_files = {}
for directory in config['search-dirs']:
    bams, indexes = find_files(directory)
    for bam in bams:
        sample = get_sample(bam.name.split('.')[0])
        if sample not in sample_files:
            sample_files[sample] = {}
        if 'bam' not in sample_files[sample]:
            sample_files[sample]['bam'] = str(bam.resolve())
    
    for index in indexes:
        sample = get_sample(index.name.split('.')[0])
        if sample in sample_files and 'bai' not in sample_files[sample]:
            sample_files[sample]['bai'] = str(index.resolve())

for _, row in files.iterrows():
    sample = row["sample_id"]
    index = row["path"].replace("bam", "bai")
    if os.path.isfile(row['path']) and os.path.isfile(index):
        if sample not in sample_files:
            sample_files[sample] = {}
        if 'bam' not in sample_files[sample]:
            sample_files[sample]['bam'] = os.path.abspath(row['path'])
            sample_files[sample]['bai'] = os.path.abspath(index)

############
# Pedigree #
############
ped = pandas.read_csv(
    config['pedigree'],
    sep="\t",
    names=["family_id", "s", "father", "mother", "sex", "is_affected"],
)

# Ensure all bam files required for the workflow are present
family_ids = set(ped.family_id)
excluded_family_ids = set()

ped_idx = ped.set_index("s", drop=False)
has_unfound_bam = 0
for s in ped['s'].unique():
    if not (s in sample_files
            and 'bam' in sample_files[s]
            and 'bai' in sample_files[s]
    ):
        logger.error("missing input files for sample {}".format(s))
        has_unfound_bam += 1
        excluded_family_ids |= set([ ped_idx.loc[s].family_id ])

family_ids = family_ids - excluded_family_ids
if has_unfound_bam > 0:
    logger.warning("excluded {} families due to missing bams ({})".format(
        has_unfound_bam, excluded_family_ids
    ))

# Ensure all samples are found in the base vcf file
# radar_vcf = pysam.VariantFile(config['cohort-vcf']['radar'])
# inova_vcf = pysam.VariantFile(config['cohort-vcf']['inova'])
cru_vcf = pysam.VariantFile(config['cohort-vcf']['cru'])
found_vcf_samples = set([
    # *radar_vcf.header.samples,
    # *inova_vcf.header.samples,
    *cru_vcf.header.samples,
])

# radar_vcf.close()
# inova_vcf.close()
cru_vcf.close()

missing_vcf_samples = set(ped_idx["s"]) - found_vcf_samples
if len(missing_vcf_samples) > 0:
    logger.error("missing vcf data for {} samples".format(len(missing_vcf_samples)))
    excluded_family_ids = set([*map(lambda s: ped_idx.loc[s].family_id, missing_vcf_samples)])
    if len(excluded_family_ids) > 0:
        logger.warning("excluding {} families due to missing vcf data ({})".format(
            len(excluded_family_ids), excluded_family_ids
        ))
        family_ids = family_ids - excluded_family_ids

#############
# Wildcards #
#############
wildcard_constraints:
    family_id="|".join(family_ids)

#############
# Functions #
#############
def get_bam(wildcards):
    trio = ped[ped.family_id == wildcards.family_id]
    samples = list(set(trio.s))
    return {
        s: sample_files[s]['bam']
        for s in samples
    }

def get_sample_bam(wildcards):
    return {
        wildcards.sample: sample_files[wildcards.sample]['bam'],
    }

def get_cohort_vcf(wildcards):
    if wildcards.family_id.startswith("R_"):
        return { 'vcf': config['cohort-vcf']['radar'] }
    elif wildcards.family_id.startswith("UTRI_"):
        return { 'vcf': config['cohort-vcf']['cru'] }
    else:
        return { 'vcf': config['cohort-vcf']['inova'] }

def get_whatshap_input(wildcards):
    return {
        'bam': 'output/{0}/{0}.bam'.format(wildcards.family_id),
        'vcf': 'output/{0}/{0}.sites.vcf.gz'.format(wildcards.family_id),
        'ped': 'output/{0}/{0}.ped'.format(wildcards.family_id),
    }

def get_unfazed_input(wildcards):
    bams = {}
    for _, row in ped[ped.family_id == wildcards.family_id].iterrows():
        # if row.father != "0" and row.mother != "0":
        bams[row.s] = "output/{0}/{0}.bam".format(row.s)
            # bams[row.s] = sample_files[row.s]['bam']
    return {
        **bams,
        'vcf': 'output/{0}/{0}.dnm_sites.vcf.gz'.format(wildcards.family_id),
        'sites': 'output/{0}/{0}.sites.vcf.gz'.format(wildcards.family_id),
        'ped': 'output/{0}/{0}.ped'.format(wildcards.family_id),
    }