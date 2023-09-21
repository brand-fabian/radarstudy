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
    raise Exception("Sample name (Output file prefix) is required.")
self_normalize = snakemake.params.get("self_normalize", False)

output_dir = os.path.dirname(snakemake.output.get("vcf"))

##
# Branch two cases:
#
#   1. Panel of Normals processing
#
if not self_normalize:
    case = snakemake.input.get("case", None)
    if case is None:
        raise Exception("Please provide a case sample.")
    normals = snakemake.input.get("normals", None)
    if normals is None:
        raise Exception("Normal input is required.")
    elif (not isinstance(normals, list)) or len(normals) < 1:
        raise Exception("Please provide an array with at least one element.")
    normals_list = " --cnv-normals-file ".join(normals)
    normals_list = " --cnv-normals-file " + normals_list

    intermediate_dir = snakemake.output.get(
        "intermediate",
        "{}-intermediate".format(output_dir)
    )
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    shell("""
        /opt/edico/bin/dragen -r {ref_dir} {extra} \
            --intermediate-results-dir {intermediate_dir} \
            --output-directory {output_dir} \
            --output-file-prefix {output_file_prefix} \
            --enable-map-align false \
            --enable-map-align-output false \
            --enable-cnv true \
            --cnv-input {case} \
            --cnv-enable-gcbias-correction false \
            {normals_list} \
            --lic-server {license_txt} {log}
    """)
else:
    ##
    #    2. Self normalization based processing [recommended for WGS]
    #
    bam = snakemake.input.get("bam", None)
    if bam is None:
        raise Exception("No bam input provided for self normalization.")
    shell("""
        /opt/edico/bin/dragen -r {ref_dir} {extra} \
            --bam-input {bam} \
            --output-directory {output_dir} \
            --output-file-prefix {output_file_prefix} \
            --enable-map-align false \
            --enable-map-align-output false \
            --enable-cnv true \
            --cnv-enable-self-normalization true \
            --lic-server {license_txt} {log}
    """)
