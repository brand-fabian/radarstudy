__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

output_vcf = snakemake.output.get("vcf", None)
if output_vcf is None:
    raise Exception("Error: Please provide an output vcf file")

extra = snakemake.params.get("extra", "")

if output_vcf.endswith('gz'):
    extra += " -O z"
elif output_vcf.endswith('bcf'):
    extra += " -O b"

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

shell("""
    bcftools norm \
        {extra} \
        -o {output_vcf} \
        {snakemake.input} {log}
    tabix {output_vcf} {append}
""")