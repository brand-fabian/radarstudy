__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

output_dir = os.path.dirname(snakemake.output[0])

variants = ""
for f in snakemake.input.variants:
    variants += " --variant {}".format(f)

shell("""
    /opt/edico/bin/dragen -r {snakemake.params.reference_dir} \
        --enable-joint-genotyping true \
        --output-file-prefix {snakemake.params.family} \
        --output-directory {output_dir} \
        --pedigree-file {snakemake.input.ped} \
        {variants} \
        --lic-server {snakemake.params.license} 2>&1 > {snakemake.log}
""")
