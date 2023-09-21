##
# Snakemake wrapper script for gatk apply base recalibration workflow
# Caution: This wrapper is using GATK v4.1
# Applies BaseRecalibrator and BQSR in one step on a single thread, because
# Spark support for GATK4 is still experimental (29.07.19)

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

recalibrator_log = snakemake.log_fmt_shell(stdout=True, stderr=True)
apply_log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

known_variants = ' '.join(map(
    lambda x: '--known-sites {}'.format(x),
    snakemake.input.known_variants
))

# Make sure we support the group ($SCRATCH_DIR workflow)
md5_required = 'md5' in snakemake.output
remote_bam = snakemake.output.bam
remote_bai = snakemake.output.bai
if os.getenv('SCRATCH_DIR'):
    remote_bam = snakemake.output.bam.replace(os.getenv('SCRATCH_DIR') + '/', '')
    remote_bai = snakemake.output.bai.replace(os.getenv('SCRATCH_DIR') + '/', '')

remote_md5 = ""
if md5_required and os.getenv('SCRATCH_DIR'):
    remote_md5 = snakemake.output.md5.replace(os.getenv('SCRATCH_DIR') + '/', '')
elif md5_required:
    remote_md5 = snakemake.output.md5

shell(r"""
# Redirect stderr to log file
exec 2> >(tee -a "{snakemake.log}")

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

gatk --java-options {jvm_args} BaseRecalibrator {extra} \
    -R {snakemake.input.ref} \
    -I {snakemake.input.bam} \
    -O "$TMPDIR/recal_table.grp" \
    {known_variants} {recalibrator_log}

gatk --java-options {jvm_args} ApplyBQSR {extra} \
    -R {snakemake.input.ref} \
    -I {snakemake.input.bam} \
    --bqsr-recal-file "$TMPDIR/recal_table.grp" \
    -O {snakemake.output.bam} \
    --create-output-bam-index \
    --create-output-bam-md5 {apply_log}

if [ $? -eq 0 ] && [ -n "$VALUE" ]
then
    touch {remote_bam} {remote_bai}
    if [ -n "{remote_md5}" ]
    then
        touch "{remote_md5}"
    fi
fi
""")
