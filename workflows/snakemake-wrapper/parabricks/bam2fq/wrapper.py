__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
import os

shell.executable('/bin/bash')

bam = snakemake.input.get('bam', None)
if bam is None:
    raise Exception("No input provided")
bam = os.path.abspath(bam)

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref = os.path.abspath(ref)

fq_out = snakemake.output.get('fastqs', None)
if fq_out is None:
    raise Exception("Please provide an output fq path")
if len(fq_out) != 2:
    raise Exception("Please specify exactly two fq output paths")

threads = snakemake.threads
extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")

# Run fastq_check tool on results
run_check = snakemake.params.get("run_check", True)
check_tool = snakemake.params.get("fastq_check", "/ceph01/homedirs/brand/Projects/radarstudy/TRIO-CRU/scripts/fastq_check/fastq_check")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = os.path.join(tempdir, str(uuid4()))
    os.makedirs(os.path.dirname(prefix), exist_ok=True)

    shell(dedent("""
        set -ex

        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbrun bam2fq {extra} --num-threads {threads} --ref {ref} \
            --in-bam {bam} \
            --out-prefix {prefix} {log}
    """))

    if run_check:
        shell(dedent("""
            set -ex
            echo "[fastq_check] Performing extra check using {check_tool}" {log_append}
            {check_tool} -fastq1 {prefix}_1.fastq.gz \
                -fastq2 {prefix}_2.fastq.gz \
                -fastqOut1 {fq_out[0]} \
                -fastqOut2 {fq_out[1]} \
                -workers {snakemake.threads} {log_append}
         """))
    else:
        for i, path in enumerate(fq_out):
            index = i + 1
            out_path = os.path.abspath(path)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            shell(dedent("""
                set -ex
                mv -vf {prefix}_{index}.fastq.gz {out_path}
            """))
