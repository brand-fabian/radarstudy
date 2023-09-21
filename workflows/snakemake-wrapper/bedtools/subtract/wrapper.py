__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

shell("""
    grep '#' {snakemake.input.a} > {snakemake.output[0]} 2>{snakemake.log}
    bedtools subtract -a {snakemake.input.a} -b {snakemake.input.b} {snakemake.params.extra} >> {snakemake.output[0]} 2>>{snakemake.log}
""")