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
output_file_prefix = snakemake.params.get("sample", None)
if output_file_prefix is None:
    raise Exception("Sample name/Output file prefix is required.")

counts = snakemake.input.get("counts", None)
if counts is None:
    raise Exception("Normalized counts input is required.")
elif (not isinstance(counts, list)) or len(counts) < 1:
    raise Exception("Please provide an array with at least one element.")
counts_list = " --cnv-input " + " --cnv-input ".join(counts)

pedigree = snakemake.input.get("pedigree", None)
if pedigree is None:
    raise Exception("Please provide a pedigree input.")


output_dir = os.path.dirname(snakemake.output.get("vcf"))

shell("""
    /opt/edico/bin/dragen -r {ref_dir} {extra} \
        --output-directory {output_dir} \
        --output-file-prefix {output_file_prefix} \
        --enable-cnv true \
        {counts_list} \
        --pedigree-file {pedigree} \
        --lic-server {license_txt} {log}
""")
