import pandas
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.logging import logger
from pathlib import Path
from collections import namedtuple
from itertools import chain
import os
import json
from meta.multiple_metrics import detect_multiple_metrics_type, get_multiple_metrics_output, METRICS_SOURCE

##
# Snakemake setup
#
report: "../report/quality-control.rst"

##
# Configuration Files
#
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##
# Setup scratch dir handling
#
work_dir = os.getenv("SCRATCH_DIR", os.getcwd())
logger.debug("Using workdir: {}".format(work_dir))

##
# Base functions
#
def find_files(base_path, file_endings=[], file_patterns=[]):
    base = Path(base_path)
    found_files = []
    for child in base.iterdir():
        # print("Found:\t{}".format(child.name))
        if child.is_file():
            file_patterns = [*file_patterns, *["*.{}".format(ending) for ending in file_endings]]
            for pattern in file_patterns:
                # print("\t{}:\t{}".format(pattern, child.match(pattern)))
                if child.match(pattern):
                    found_files.append(child)
                    break
        if child.is_dir():
            found_files = [
                *found_files,
                *find_files(str(child), file_endings=file_endings, file_patterns=file_patterns)
            ]
    return found_files

def get_cohort_vcf(base = config.get("cohort_vcf", None)):
    """Get the path to the cohort vcf.
    
    Returns
    -------
    bool
        True, if all vcfs found in the input dir should be merged
        to create the cohort vcf, False if it is derived from the
        config parameter "cohort_vcf".
    str
        Path to the output file.
    """
    if base is None or not os.path.isfile(base):
        return True, "output/{project}/{project}.vcf.gz"
    else:
        if base.endswith("bcf"):
            # Convert bcf to .vcf.gz
            return False, "output/{project}/{project}.vcf.gz"
        else:
            return False, os.path.abspath(base)

##
# Input definition
#
endings = config.get("file-endings", [])
patterns = config.get("search-patterns", [])
files = []
for path in config["input-dirs"]:
    files = [
        *files,
        *find_files(path, file_endings=endings, file_patterns=patterns)
    ]

default_metrics_source = config.get("default_metrics_source", "picard")
if default_metrics_source == "picard":
    default_metrics_source = METRICS_SOURCE.PICARD
elif default_metrics_source == "parabricks":
    default_metrics_source = METRICS_SOURCE.PARABRICKS
else:
    logger.error("Invalid metrics source {}. It must be one of {{picard, parabricks}}".format(default_metrics_source))
metrics_source = detect_multiple_metrics_type(map(str, files), default=default_metrics_source)
logger.info("Detected multiple metrics source: {}".format(metrics_source))

samples = dict()
for f in files:
    sample_name = f.name.split(".")[0]
    if "fastq" in f.name or "fq" in f.name:
        sample_name = f.name.split("_")[0]
    if not sample_name in samples.keys():
        samples[sample_name] = []
    if not any(map(lambda path: f.name == path.name, samples[sample_name])):
        samples[sample_name].append(f)

wildcard_constraints:
    sample="|".join(samples.keys()),
    project=config["project"]

def multiple_metrics_output():
    return get_multiple_metrics_output(metrics_source)

##
# Helper variables
#
alignment_stats = [
    "samtools-stats.txt",
    "flagstat.txt",
    "idxstats.txt",
    "mosdepth.global.dist.txt",
    "mosdepth.region.dist.txt",
    "mosdepth.summary.txt",
    "wgs_metrics.txt",
    "selfSM"
] + list(map(lambda x: ".".join(x.split(".")[1:]), multiple_metrics_output().values()))

variant_stats = [
    "bcftools-stats.txt"
]

cohort_stats = [
    "bcftools-stats.txt",
    "relatedness2"
]

dragen_stats = [
    "vc_metrics.csv",
    "ploidy_estimation_metrics.csv",
    "wgs_contig_mean_cov.csv",
    "wgs_coverage_metrics.csv",
    "wgs_fine_hist.csv",
    "fragment_length_hist.csv",
    "mapping_metrics.csv",
]

##
# Helper functions
#
DefaultWildcard = namedtuple('DefaultWildcard', [ 'sample' ])

def get_fastq(wildcards):
    return list(map(str, filter(lambda f: any(str(f).endswith(e) for e in [
        "fastq", "fq", "fastq.gz", "fq.gz"
    ]), samples[wildcards.sample])))

def get_bams(wildcards):
    bams = list(
        map(
            str,
            filter(lambda f: "bam" in f.name and not "bai" in f.name, samples[wildcards.sample])
        )
    )
    bais = list(range(len(bams)))
    for idx, bam in enumerate(bams):
        bai = "{}.bai".format(bam)
        alternative = ".".join([*bai.split(".")[:-2], "bai"])
        if not Path(bai).exists() and not Path(alternative).exists():
            raise Exception(
                "Could not find bam index {} or alternative {}.".format(
                    bai,
                    alternative
                )
            )
        elif Path(bai).exists():
            bais[idx] = bai
        elif Path(alternative).exists():
            bais[idx] = alternative
    return {
        "bam": bams,
        "bai": bais
    }

def get_bam(wildcards):
    all_bams = get_bams(wildcards)
    if len(all_bams["bam"]) and len(all_bams["bai"]):
        return {
            "bam": all_bams["bam"][0],
            "bai": all_bams["bai"][0]
        }
    else:
        return all_bams

def get_variants(wildcards):
    """Find and return all .vcf (variant call format) files found in the
       initial search directories."""
    return [str(x) for x in filter(lambda name: name.match("*vcf*"), files)]

def get_variant(wildcards):
    """Find and return all .vcf files found for the given sample."""
    return list(map(str, filter(lambda f: f.match("*vcf*"), samples[wildcards.sample])))

def get_ngs_chew_fingerprints(wildcards):
    vcf = []
    npz = []
    for sample in samples.keys():
        alignment = get_bams(DefaultWildcard(sample=sample))
        if len(alignment["bam"]) and len(alignment["bai"]):
            for f_name in alignment["bam"]:
                vcf.append("output/{sample}/{sample}.fingerprint.vcf.gz".format(
                    sample=sample
                ))
                npz.append("output/{sample}/{sample}.npz".format(sample=sample))
    return {
        "vcf": vcf,
        "npz": npz
    }

def get_multiqc_data(wildcards):
    files = []
    for sample in samples.keys():
        for idx, f_name in enumerate(get_fastq(DefaultWildcard(sample=sample))):
            yield "output/{sample}/{sample}_{idx}_fastqc.zip".format(
                sample=sample,
                idx=idx + 1
            )
        alignment = get_bams(DefaultWildcard(sample=sample))
        if len(alignment["bam"]) and len(alignment["bai"]):
            for idx, bam in enumerate(alignment["bam"]):
                for ending in alignment_stats:
                    yield "output/{sample}/{sample}.{ending}".format(
                        sample=sample,
                        ending=ending.format(sample=sample),
                    )
        for idx, variant in enumerate(get_variant(DefaultWildcard(sample=sample))):
            if len(variant):
                for ending in variant_stats:
                    yield "output/{sample}/{sample}.{ending}".format(
                        sample=sample,
                        ending=ending.format(sample=sample),
                    )

    if config.get("cohort_vcf", None) is not None:
        for ending in cohort_stats:
            yield "output/{project}/{project}.{ending}".format(
                project=config["project"],
                ending=ending,
            )

    if config.get("pedigree", None) is not None:
        ##
        # Execute peddy rules, if pedigree is supplied as input
        yield "output/{project}/{project}.peddy.ped".format(
            project=config["project"]
        )
        yield "output/{project}/{project}.ped_check.csv".format(
            project=config["project"]
        )
        yield "output/{project}/{project}.sex_check.csv".format(
            project=config["project"]
        )
        yield "output/{project}/{project}.background_pca.json".format(
            project=config["project"]
        )
        yield "output/{project}/{project}.het_check.csv".format(
            project=config["project"]
        )

    if config.get("dragen_qc", False):
        # Dragen QC is expected to be part of the input
        for f in files:
            if any(f.endswith(e) for e in dragen_stats):
                yield f
