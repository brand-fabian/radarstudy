##
# Snakemake wrapper script for bwa mem
# Adapted from https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bwa/mem.html
# Original authors: Johannes Köster, Julian de Ruiter

__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
sort_order = snakemake.params.get("sort_order", "coordinate")
sort_extra = snakemake.params.get("sort_extra", "")
sort_jvm_args = snakemake.params.get("sort_jvm_args", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if not isinstance(snakemake.input.fq, str) and len(snakemake.input.fq) % 2 != 0:
    raise ValueError("input must have 1 (single-end) or "
                     "2 (paired-end) elements")

if sort_order not in {"coordinate", "queryname"}:
    raise ValueError("Unexpected value for sort_order ({})".format(sort_order))

# Make sure we support the group ($SCRATCH_DIR workflow)

md5_required = 'md5' in snakemake.output
remote_bam = snakemake.output.bam
remote_bai = snakemake.output.bai
if os.getenv('SCRATCH_DIR'):
    remote_bam = snakemake.output.bam.replace(os.getenv('SCRATCH_DIR') + '/', '')
    if not os.path.exists(os.path.dirname(remote_bam)):
        os.makedirs(os.path.dirname(remote_bam))
    remote_bai = snakemake.output.bai.replace(os.getenv('SCRATCH_DIR') + '/', '')
    if not os.path.exists(os.path.dirname(remote_bai)):
        os.makedirs(os.path.dirname(remote_bai))

remote_md5 = ""
if md5_required and os.getenv('SCRATCH_DIR'):
    remote_md5 = snakemake.output.md5.replace(os.getenv('SCRATCH_DIR') + '/', '')
elif md5_required:
    remote_md5 = snakemake.output.md5

##
# Wrapper shell code to execute bwa mem and picard SortSam
# To support execution on local drive, we touch the file path
# of the remote machine if execution was successful.
shell(r"""
# Redirect stderr to log file
exec 2> >(tee -a "{snakemake.log}")

exec_align_sort () {{
    if [ -z "$1" ]; then
        echo "exec_align_sort expects one argument: OUTPUT_FILE"
    fi

    (
        bwa mem -t {snakemake.threads} \
            {extra} \
            {snakemake.params.index} \
            {snakemake.input.fq} | picard {sort_jvm_args} \
                SortSam {sort_extra} \
                TMP_DIR=$TMPDIR \
                INPUT=/dev/stdin \
                OUTPUT=$1 \
                SORT_ORDER={sort_order} \
                CREATE_INDEX=true \
                CREATE_MD5_FILE=true
    ) {log}
}}

export TARGET_BAM="{snakemake.output.bam}"

VALUE=${{SCRATCH_DIR:-}}
if [ -z "$VALUE" ]; then
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    echo "Running alignment and sorting on shared drive with tmp dir $TMPDIR"
    exec_align_sort $TARGET_BAM
else
    export TMPDIR=$SCRATCH_DIR/bwa
    mkdir -p "$TMPDIR"
    mkdir -p $(dirname "$TARGET_BAM")

    echo "Running alignment and sorting on local scratch space $SCRATCH_DIR"
    exec_align_sort $TARGET_BAM

    # Touch files on remote, if possible
    if [ $? -eq 0 ]; then
        touch "{remote_bam}" "{remote_bai}"
        if [ -n "{remote_md5}" ]
        then
            touch "{remote_md5}"
        fi
    fi
fi
""")
