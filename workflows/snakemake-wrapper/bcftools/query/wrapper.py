__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

output_file = snakemake.output.get("file", None)
if output_file is None:
    raise Exception("Error: Please provide an output file file")

format_str = snakemake.params.get('format', None)
if format_str is None:
    raise Exception("Error: Please provide a format string to query.")

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("""
    bcftools query \
        {extra} \
        -f '{format_str}' \
        {snakemake.input} {log} > {output_file}
""")