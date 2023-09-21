##
# Snakemake wrapper script for peddy
# Runs the tool in a temp dir and only moves desired output files
# to target paths (similar to the behaviour of a shadow rule)
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
import uuid
from tempfile import TemporaryDirectory

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
threads = snakemake.threads if hasattr(snakemake, 'threads') else 1

known_suffix = {
    "ped": "peddy.ped",
    "het_check": "het_check.csv",
    "ped_check": "ped_check.csv",
    "sex_check": "sex_check.csv",
    "background_pca": "background_pca.json",
    "html": "html"
}

prefix = snakemake.params.get("prefix", str(uuid.uuid4()))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory() as tempdir:
    output_prefix = os.path.join(tempdir, prefix)
    shell(r"""
        echo "Target prefix: {output_prefix}"
        peddy -p {snakemake.threads} {extra} --plot \
            --prefix {output_prefix} \
            {snakemake.input.vcf} \
            {snakemake.input.ped} {log}
        ls -l {tempdir}
    """)

    ##
    # Move desired output files to target location
    for key, value in known_suffix.items():
        if key in snakemake.output.keys():
            source_path = output_prefix + "." + value
            target_path = snakemake.output.get(key)
            if target_path != source_path:
                shell("mv {source_path} {target_path}")