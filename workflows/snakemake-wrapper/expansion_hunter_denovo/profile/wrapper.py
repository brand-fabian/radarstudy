__author__ = "Fabian Brand"
__email__ = "brand@imbie.uni-bonn.de"

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import uuid
import os

shell.executable("/bin/bash")

CATALOG_HG19 = "https://github.com/Illumina/RepeatCatalogs/raw/master/hg19/variant_catalog.json"
CATALOG_HG38 = "https://github.com/Illumina/RepeatCatalogs/raw/master/hg38/variant_catalog.json"

output_files = {
    "profile": ".str_profile.json",
    "locus": ".locus.tsv",
    "motif": ".motif.tsv",
    "reads": ".reads.tsv",
}

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("No input reads path set.")
if not os.path.isfile(bam):
    raise Exception("Missing input bam {}.".format(bam))
bam_str = "--reads " + os.path.abspath(bam)

ref = snakemake.params.get("reference", None)
if ref is None:
    raise Exception("No reference genome path supplied.")
if not os.path.isfile(ref):
    raise Exception("Reference genome {} does not exist.".format(ref))
ref_str = "--reference " + os.path.abspath(ref)

extra = snakemake.params.get("extra", "")
min_anchor_mapq = snakemake.params.get("min_anchor_mapq", 50)
max_irr_mapq = snakemake.params.get("max_irr_mapq", 40)
threads = snakemake.threads or 1
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
log_f = snakemake.log[0]

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = str(uuid.uuid4())
    shell(dedent("""
    set -ex

    echo > {log_f}

    {{
        ExpansionHunterDenovo profile \
            {extra} \
            {bam_str} \
            {ref_str} \
            --output-prefix {tempdir}/{prefix} \
            --min-anchor-mapq {min_anchor_mapq} \
            --max-irr-mapq {max_irr_mapq}
    }} {log}
    """))

    for key, target in snakemake.output.items():
        if key in output_files:
            source = output_files[key]
            shell(dedent("""
                set -ex
                {{
                    mv -v {tempdir}/{prefix}{source} {target}
                }} {log}
            """))
