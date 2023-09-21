__author__ = "Fabian Brand"
__email__ = "brand@imbie.uni-bonn.de"

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import uuid
import os

shell.executable("/bin/bash")

manifest = snakemake.input.get("manifest", None)
if manifest is None or not os.path.isfile(manifest):
    raise Exception("No or missing manifest file supplied.")
manifest_str = "--manifest " + os.path.abspath(manifest)

multisample_profile = snakemake.input.get("multisample_profile", None)
if multisample_profile is None or not os.path.isfile(multisample_profile):
    raise Exception("No or missing multisample profile supplied.")
multisample_profile_str = "--multisample-profile " + os.path.abspath(multisample_profile)

output_tsv = snakemake.output.get("tsv", None)
if output_tsv is None:
    raise Exception("Output filename missing.")
output_str = "--output " + os.path.abspath(output_tsv)

extra = snakemake.params.get("extra", "")
analysis_type = snakemake.params.get("analysis_type", "locus")
threads = snakemake.threads or 1
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
log_f = snakemake.log[0]

assert analysis_type in ["locus", "motif"]

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = str(uuid.uuid4())
    shell(dedent("""
    set -ex

    echo > {log_f}

    {{
        outlier.py {analysis_type} \
            {extra} \
            {manifest_str} \
            {multisample_profile_str} \
            {output_str}
    }} {log}
    """))

