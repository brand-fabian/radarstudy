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
    "json": ".json",
    "vcf": ".vcf",
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

catalog = snakemake.params.get("catalog", CATALOG_HG19)
extra = snakemake.params.get("extra", "")
analysis_mode = snakemake.params.get("analysis_mode", "streaming")
sex = snakemake.params.get("sex", "female")
module_name = snakemake.params.get("module_name", "bio/ExpansionHunter")
threads = snakemake.threads or 1
log_reset = snakemake.log_fmt_shell(stdout=True, stderr=True)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = str(uuid.uuid4())
    shell("echo {log_reset}")
    if catalog.startswith("http") or catalog.startswith("ftp"):
        shell(dedent("""
            set -ex
            {{
                curl -L {catalog} > {tempdir}/{prefix}.catalog.json
            }} {log}
        """))
    else:
        if not os.path.isfile(catalog):
            raise Exception("Catalog {} is not downloadable and does not exist as file.".format(catalog))
        shell(dedent("""
            set -ex
            {{
                ln -s {catalog} {tempdir}/{prefix}.catalog.json
            }} {log}
        """))
    
    shell(dedent("""
        set -ex
        # Load the expansion hunter module if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        {{
            ExpansionHunter {extra} \
                --threads {threads} \
                {bam_str} \
                --variant-catalog {tempdir}/{prefix}.catalog.json \
                {ref_str} \
                --analysis-mode {analysis_mode} \
                --sex {sex} \
                --output-prefix {tempdir}/{prefix}
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
