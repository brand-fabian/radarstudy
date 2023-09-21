##
# Snakemake wrapper script for gatk haplotype caller
# Caution: This wrapper is using GATK v4.1
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

if not jvm_args == "" and jvm_args is not None:
    jvm_args = "--java-options " + jvm_args

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

dbsnp = snakemake.input.get("dbsnp", "")
if dbsnp:
    dbsnp = "--dbsnp " + dbsnp

interval = snakemake.input.get("interval", "")
if not isinstance(interval, str):
    if len(interval) > 1:
        interval = " --intervals ".join(interval)
        interval = "--intervals " + interval
        if "--interval-set-rule" not in extra and "-isr" not in extra:
            interval += " --interval-set-rule INTERSECTION"
    else:
        raise Exception("Please specify intervals either as string, or array of strings.")
elif not interval == "":
    interval = " --intervals " + interval

# Make sure we support the group ($SCRATCH_DIR workflow)
remote_gvcf = snakemake.output.gvcf
remote_idx = snakemake.output.idx
if os.getenv('SCRATCH_DIR'):
    remote_gvcf = snakemake.output.gvcf.replace(os.getenv('SCRATCH_DIR'), '')
    remote_idx = snakemake.output.idx.replace(os.getenv('SCRATCH_DIR'), '')

shell(r"""
# Redirect stderr to log file
VALUE=${{SCRATCH_DIR:-}}

if [ -z "$VALUE" ]
then
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT
    echo "Using temp dir $TMPDIR."
else
    export TMPDIR=$SCRATCH_DIR/gatk
    mkdir -p "$TMPDIR"
    echo "Using temp dir $TMPDIR."
fi

gatk {jvm_args} HaplotypeCaller \
    {extra} \
    -R {snakemake.input.ref} \
    -I {snakemake.input.bam} \
    -O {snakemake.output.gvcf} \
    {interval} \
    -ERC GVCF \
    --create-output-variant-index \
    --tmp-dir $TMPDIR \
    {dbsnp} {log}

if [ $? -eq 0 ] && [ -n "$VALUE" ]
then
    touch {remote_gvcf} {remote_idx}
fi
""")
