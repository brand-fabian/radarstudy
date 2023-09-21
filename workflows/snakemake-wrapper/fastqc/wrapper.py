##
# Snakemake wrapper for fastqc
# Built upon: https://bitbucket.org/snakemake/snakemake-wrappers/src/master/bio/fastqc/wrapper.py

__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from tempfile import TemporaryDirectory, gettempdir

def basename_without_ext(file_path):
    split_ind = 2 if file_path.endswith('.gz') else 1
    return ".".join(os.path.basename(file_path).split(".")[:-split_ind])

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    print("Using temp dir {}".format(tempdir))
    shell(r"""
    fastqc {extra} \
        --quiet \
        --outdir "{tempdir}" \
        {snakemake.input} \
        {log}
    """)

    ##
    # Move files to target location
    if (len(snakemake.input) == len(snakemake.output.html)
            and len(snakemake.input) == len(snakemake.output.zip)):
        for idx, path in enumerate(snakemake.input):
            base = basename_without_ext(path)
            zip_path = os.path.join(tempdir, base + "_fastqc.zip")
            html_path = os.path.join(tempdir, base + "_fastqc.html")
            zip_target = snakemake.output.zip[idx]
            html_target = snakemake.output.html[idx]

            if zip_path != zip_target:
                shell("mv {zip_path} {zip_target}")
            if html_path != html_target:
                shell("mv {html_path} {html_target}")
    else:
        raise Exception("More input files than output files for fastqc.")

    ##
    # Mark output, if using scratch dir
    if os.getenv('SCRATCH_DIR'):
        for path in [*snakemake.output.html, *snakemake.output.zip]:
            remote = path.replace(os.getenv('SCRATCH_DIR'), '')
            shell("touch {remote}")