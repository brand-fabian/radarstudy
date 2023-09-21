__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("stats", "")

shell(r"""
    bcftools stats {extra} \
        {snakemake.input} > {snakemake.output} {log}
""")