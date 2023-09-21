__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import pandas
import os

FLANKING_BASES = 1000

vcf_in = snakemake.input.get("vcf", None)
if vcf_in is None:
    raise Exception("no vcf input provided")
vcf_in = os.path.abspath(vcf_in)

vcf_out = snakemake.output.get("vcf", None)
if vcf_out is None:
    raise Exception("no vcf output provided")

ped = snakemake.params.get("pedigree", None)
if ped is None:
    raise Exception("no pedigree file provided")
ped = os.path.abspath(ped)

sites = snakemake.params.get("sites", None)
if sites is None:
    raise Exception("no sites table provided")
sites = os.path.abspath(sites)

family_id = snakemake.params.get("family_id", None)
if family_id is None:
    raise Exception("no family_id provided")

flanking_bases = snakemake.params.get("flanking_bases", FLANKING_BASES)

ped_df = pandas.read_csv(
    ped,
    sep="\t",
    names=["family_id", "s", "father", "mother", "sex", "is_affected"]
)
trio = ped_df[ped_df.family_id == family_id]
samples = set(trio.s)
dnm_df = pandas.read_csv(sites, sep="\t")
dnm_df = dnm_df[dnm_df.s.isin(samples)]

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    regions_fname = os.path.join(tempdir, "dnm.regions")
    with open(regions_fname, 'w') as regions_f:
        for _, row in dnm_df.iterrows():
            contig, pos = row.locus.split(':')
            regions_f.write("\t".join([
                contig,
                str(max(0, int(pos) - flanking_bases)),
                str(int(pos) + flanking_bases)
            ]) + "\n")

    sample_list = ",".join(samples)

    if len(dnm_df) > 0:
        regions_cli = "-R {}".format(regions_fname)
    else:
        regions_cli = ""

    shell(dedent("""
        set -x
        cp {regions_fname} output/{family_id}/{family_id}.regions
        bcftools view -a -U -c 1 {regions_cli} \
            --force-samples -s {sample_list} \
            {vcf_in} | bgzip > {vcf_out}
        tabix {vcf_out}
    """))
