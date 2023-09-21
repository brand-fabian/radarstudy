##
# Snakemake wrapper script for picard MergeVcfs
#
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

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

input_files = ""
for in_file in snakemake.input.vcfs:
    input_files += ("I=" + in_file + " ")

##
# Make sure we support the group ($SCRATCH_DIR workflow)
remote_vcf = snakemake.output.vcf
remote_idx = snakemake.output.idx
if os.getenv('SCRATCH_DIR'):
    remote_vcf = snakemake.output.vcf.replace(os.getenv('SCRATCH_DIR'), '')
    remote_idx = snakemake.output.idx.replace(os.getenv('SCRATCH_DIR'), '')

md5_required = 'md5' in snakemake.output
remote_md5 = ""
if md5_required and os.getenv('SCRATCH_DIR'):
    remote_md5 = snakemake.output.md5.replace(os.getenv('SCRATCH_DIR') + '/', '')
elif md5_required:
    remote_md5 = snakemake.output.md5

shell(r"""

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

picard {jvm_args} MergeVcfs {extra} \
    {input_files} \
    OUTPUT={snakemake.output.vcf} \
    CREATE_INDEX=true \
    TMP_DIR=$TMPDIR \
    CREATE_MD5_FILE=true {log}

if [ $? -eq 0 ] && [ -n "$VALUE" ]
then
    touch {remote_vcf} {remote_idx} 
    if [ -n "{remote_md5}" ]
    then
        touch {remote_md5}
    fi
fi
""")
