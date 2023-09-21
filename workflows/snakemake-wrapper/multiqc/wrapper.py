##
# Snakemake wrapper script for gatk SelectVariants
# Caution: This wrapper is using GATK v4.1
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
import uuid
from tempfile import TemporaryDirectory

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")

input_dirs = set(os.path.dirname(f) for f in snakemake.input)
output_dir = os.path.dirname(snakemake.output[0])
output_name = os.path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(r"""
    multiqc {extra} \
        --force \
        -o {output_dir} \
        -n {output_name} \
        {input_dirs} {log}
""")