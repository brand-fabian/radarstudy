__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

license_txt = snakemake.params.get("license", None)
if license_txt is None:
    raise Exception("License is provided for this job.")
ref_dir = snakemake.params.get("reference_dir", None)
if ref_dir is None:
    raise Exception("Please provide the path to a cnv calling enabled dragen hash table.")

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("Bam input is required.")
output_file_prefix = snakemake.params.get("sample", None)
if output_file_prefix is None:
    raise Exception("Sample name is required.")
output_dir = os.path.dirname(snakemake.output.get("normal"))
intermediate_dir = snakemake.output.get(
    "intermediate",
    "{}-intermediate".format(output_dir)
)
if not os.path.exists(intermediate_dir):
    os.makedirs(intermediate_dir)

shell("""
    /opt/edico/bin/dragen -r {ref_dir} {extra} \
        -b {bam} \
        --intermediate-results-dir {intermediate_dir} \
        --output-directory {output_dir} \
        --output-file-prefix {output_file_prefix} \
        --enable-map-align false \
        --enable-map-align-output false \
        --enable-cnv true \
        --lic-server {license_txt} {log}
""")
