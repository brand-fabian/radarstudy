__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

shell("""
    python filter/main.py {snakemake.params.extra} --output-vcf {snakemake.output[0]} {snakemake.input} {snakemake.params.sample} 2>{snakemake.log}
""")