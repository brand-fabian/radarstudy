__author__ = 'Sugi Sivalingam'
__email__ = 's.sivalingam@uni-bonn.de'


# check file existence
# -M as defined parameter (default 1e-8)
# options -Oz (check file ending - gz?!)
# seperate snakemake input (vcf / gnomad ref)

from snakemake.shell import shell


shell.executable('/bin/bash')
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

shell(r"""
    bcftools roh {extra} --AF-file {snakemake.input.gnomad} --threads {snakemake.threads} {snakemake.input.vcf} -o {snakemake.output} {log}
""")
