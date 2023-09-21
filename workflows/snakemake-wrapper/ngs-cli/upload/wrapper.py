##
# Snakemake wrapper script for ngs-cli upload workflow
# Uploads all input files to the ngs-platform
#
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

project = snakemake.params.get("project", "")
sample = snakemake.params.get("sample", None)

if sample is None:
    raise Exception("Please provide a sample id to upload files.")

num_inputs = len(snakemake.input)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
shell(r"""
# Redirect stderr to log file
set +e
exec 2> >(tee -a "{snakemake.log}")

if [ -n "{project}" ]
then
    OLD_PRJ=$(ngs-cli config ls | grep 'Project' | awk -F ' ' '{{print $2}}')
    ngs-cli config set-project "{project}"
    if [ -n "$OLD_PRJ" ]
    then
        trap "ngs-cli config set-project \"$OLD_PRJ\"" EXIT
    fi
fi

SID=$(ngs-cli sample find "{sample}")
if [[ "$SID" =~ "^Could not.*" ]]; then
    echo "Could not find sample {snakemake.params.sample}"
    exit 1
fi

ngs-cli files upload $SID {snakemake.input:q} {log} | grep -qi "error"

if [ $? = 0 ]; then
    echo "Error uploading {snakemake.input}"
    exit 1
else
    touch {snakemake.output:q}
fi
""")
