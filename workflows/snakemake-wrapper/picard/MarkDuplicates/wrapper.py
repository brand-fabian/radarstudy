##
# Snakemake wrapper script for picard mark duplicates taking into account
# local requirements, such as local scratch dirs on the cluster.

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

##
# Make sure we support the group ($SCRATCH_DIR workflow)
# Touch one remote file at the same path cwd + output 
# instead of $SCRATCH_DIR + output
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

picard {jvm_args} MarkDuplicates {extra} \
    INPUT={snakemake.input.bam} \
    OUTPUT={snakemake.output.bam} \
    METRICS_FILE={snakemake.output.metrics} \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true {log}

if [ $? -eq 0 ] && [ -n "$VALUE" ]
then
    touch {remote_bam} {remote_bai}
    if [ -n "{remote_md5}" ]
    then
        touch "{remote_md5}"
    fi
fi
""")
