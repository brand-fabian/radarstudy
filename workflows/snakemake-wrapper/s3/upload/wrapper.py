__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

sample_name = snakemake.params.get("sample", None)
if sample_name is None:
    raise Exception("Sample name is required to upload sample data.")

input_dir = os.path.dirname(snakemake.input[0])
output_dir = os.path.dirname(snakemake.output[0])
shell("""
    aws s3 sync {input_dir} {output_dir} --exclude \"*\" --include \"{sample_name}*\" {log}
""")
