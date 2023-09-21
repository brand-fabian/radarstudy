import pandas
from rules.meta.file import ensure_bams
from snakemake.utils import validate, min_version

from meta import ensure_fastqs, get_multiple_metrics_output, METRICS_SOURCE

##
# Setup snakemake and read configuration files
#
min_version("6.1.0")
report: "../report/workflow.rst"

configfile: "config.yml"
validate(config, schema="../schemas/config.schema.yml")

samples = pandas.read_csv(
    config["sample-sheet"],
    sep="\t",
    header=0,
).set_index("sample_id", drop=False)
validate(samples, schema="../schemas/samples.schema.yml")

files = pandas.read_csv(
    config["files-table"],
    sep="\t",
    header=0,
).set_index("sample_id", drop=False)
validate(files, schema="../schemas/files.schema.yml")


fastqs = ensure_fastqs(files, samples)
bams = ensure_bams(files, samples)


wildcard_constraints:
    sample="|".join(samples.index)

##
# Resolve file paths
#
def get_bam(sample):
    if sample in bams and len(bams[sample]) > 0:
        return bams[sample][0]
    else:
        return None


def get_fq2bam_input(wildcards):
    assert wildcards.sample in fastqs and len(fastqs[wildcards.sample]) > 0
    return {
        "fq": sorted(fastqs[wildcards.sample])
    }


def get_alignment(wildcards):
    input_bams = get_bam(wildcards.sample)
    if input_bams is not None:
        bam, bai = input_bams
        return {
            "bam": bam,
            "bai": bai,
        }
    else:
        return {
            "bam": "output/{sample}/{sample}.bam".format(sample=wildcards.sample),
            "bai": "output/{sample}/{sample}.bam.bai".format(sample=wildcards.sample),
        }
