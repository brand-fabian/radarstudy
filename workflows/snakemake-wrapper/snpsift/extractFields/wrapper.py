__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

input_vcf = snakemake.input.get('vcf', None)
if input_vcf is None:
    raise Exception("Please provide an input vcf.")
fields = set(snakemake.params.get('fields', None))
if fields is None:
    raise Exception("Please provide a set of fields to extract.")

output = snakemake.output
if output is None:
    raise Exception("Please provide an output vcf path")

if input_vcf.endswith('gz'):
    input_str = f"zcat {input_vcf}"
else:
    input_str = f"cat {input_vcf}"

fields_str = ' '.join(fields)

output_str = f"> {output}"

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("""
    ({input_str} | SnpSift extractFields {extra} - {fields_str} {output_str}) {log}
""")