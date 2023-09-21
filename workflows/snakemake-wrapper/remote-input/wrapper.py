##
# Snakemake wrapper script for ngs-cli download workflow

__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
remote = snakemake.params.get("remote", None)
project = snakemake.params.get("project", None)
bucket = snakemake.params.get("bucket", None)
sample = snakemake.params.get("sample", None)

if remote is None:
    raise Exception("Please provide a remote type to get input files.")

if len(snakemake.input) == len(snakemake.output):
    raise Exception("Invalid Setup: Different number of input and output files.")

if remote == "s3" and bucket is None:
    raise Exception("Error: Remote is S3 and no bucket provided.")

if remote == "ngs-cli" and project is None:
    raise Exception("Error: Remote is ngs-cli and no project provided.")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

if remote == "s3":
    # Copy files from remote
    for idx, in_file in enumerate(snakemake.input.items()):
        output_file = snakemake.output.get(in_file[0], None)
            if output_file is None:
                output_file = snakemake.output[idx]
        shell("aws s3 cp s3://{}/{} {} {}".format(
            bucket,
            in_file[1],
            output_file,
            log if idx == 0 else log_append
        ))
elif remote == "ngs-cli":
    # Copy files from ngs-cli to target path
    for idx, in_file in enumerate(snakemake.input.items()):
        output_file = snakemake.output.get(in_file[0], None)
            if output_file is None:
                output_file = snakemake.output[idx]
        shell(r"""
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
FID=$(ngs-cli files ls --like "$(basename {in_file[1]})" $SID | tail -1 | awk -F ' ' '{{print $1}}')
ngs-cli files download --target "{snakemake.output[0]}" $SID $FID {log_append}
        """)
else:
    # Link files to output for local files
    for idx, in_file in enumerate(snakemake.input.items()):
        output_file = snakemake.output.get(in_file[0], None)
            if output_file is None:
                output_file = snakemake.output[idx]
        shell("ln -s {} {} {}".format(
            os.path.abspath(in_file[1]),
            os.path.abspath(output_file),
            log if idx == 0 else log_append
        )
