##
# Snakemake wrapper script for ngs-cli download workflow
# Limitation: Max. of 2 files is downloaded
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
project = snakemake.params.get("project", "")
sample = snakemake.params.get("sample", None)

if sample is None:
    raise Exception("Please provide a sample to download files.")

num_outputs = len(snakemake.output)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(r"""
# Redirect stderr to log file
exec 2> >(tee -a "{snakemake.log}")

TEST=${{SCRATCH_DIR:-}}
if [ -z "$TEST" ]
then
    echo "Downloading to remote locations..."
else
    echo "Downloading to local scratch $SCRATCH_DIR"
fi

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
    echo "Could not find sample {sample}"
    exit 1
fi

if [ ! -d "$(dirname {snakemake.output[0]})" ]
then
    mkdir -p "$(dirname {snakemake.output[0]})"
fi

FID=$(ngs-cli files ls --like "$(basename {snakemake.input[0]})" $SID | tail -1 | awk -F ' ' '{{print $1}}')
ngs-cli files download --target "{snakemake.output[0]}" $SID $FID {log}
if [ -n "$TEST" ]
then
    TARGET={snakemake.output[0]}
    touch ${{TARGET##${{SCRATCH_DIR}}\/}}
fi

if [ {num_outputs} -gt 1 ]
then
    if [ ! -d "$(dirname {snakemake.output[0]})" ]
    then
        mkdir -p "$(dirname {snakemake.output[0]})"
    fi
    FID2=$(ngs-cli files ls --like "$(basename {snakemake.input[1]})" $SID | tail -1 | awk -F ' ' '{{print $1}}')
    ngs-cli files download --target "{snakemake.output[1]}" $SID $FID2 {log}
    if [ -n "$TEST" ]
    then
        TARGET={snakemake.output[1]}
        touch ${{TARGET##${{SCRATCH_DIR}}\/}}
    fi
fi
""")
