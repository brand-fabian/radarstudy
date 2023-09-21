__author__ = "Fabian Brand"
__email__ = "brand@imbie.uni-bonn.de"

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import uuid
import os

shell.executable("/bin/bash")

output_files = {
    "profile": ".multisample_profile.json",
    "manifest": ".manifest.tsv",
}

profiles = snakemake.input.get("profiles", [])
if len(profiles) == 0:
    raise Exception("Please supply at least one input profile.")
if any(not os.path.isfile(p) for p in profiles):
    raise Exception("Some input profiles are missing. Please make sure to provide the correct paths.")
profiles = list(map(os.path.abspath, profiles))

ref = snakemake.params.get("reference", None)
if ref is None:
    raise Exception("No reference genome path supplied.")
if not os.path.isfile(ref):
    raise Exception("Reference genome {} does not exist.".format(ref))
ref_str = "--reference " + os.path.abspath(ref)

extra = snakemake.params.get("extra", "")
sample_names = snakemake.params.get("sample_names", [])
case_control_status = snakemake.params.get("case_control", [])
threads = snakemake.threads or 1
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
log_f = snakemake.log[0]

assert len(sample_names) == len(profiles)
assert len(case_control_status) == len(profiles)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = str(uuid.uuid4())
    with open("{}/{}.manifest.tsv".format(tempdir, prefix), "w") as f:
        for row in zip(sample_names, case_control_status, profiles):
            f.write("\t".join(row) + "\n")

    shell(dedent("""
    set -ex

    echo > {log_f}

    {{
        ExpansionHunterDenovo merge \
            {extra} \
            --manifest {tempdir}/{prefix}.manifest.tsv \
            {ref_str} \
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
