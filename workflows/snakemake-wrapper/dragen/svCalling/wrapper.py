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
    raise Exception("Please provide the path to a sv calling enabled dragen hash table.")
fasta_ref = snakemake.params.get("fasta_reference", None)
if fasta_ref is None:
    fasta = ''
else:
    fasta = '--sv-reference {fasta_ref}'.format(fasta_ref=fasta_refs)

bams = snakemake.input.get("bams", None)
if bams is None:
    raise Exception("Bam input is required.")
elif (not isinstance(bams, list)) or len(bams) < 3:
    raise Exception("Please provide an array with at least three elements.")
bams_list = " --bam-input " + " --bam-input ".join(bams)

bais = snakemake.input.get("bais", None)
##
# Check that bam inputs match bais
for bam in bams:
    found = False
    for bai in bais:
        if bai.startswith(bam.split(".")[0]):
            found = True
            break
    if not found:
        raise Exception("No matching .bai index for {}".format(bam))

pedigree = snakemake.input.get("pedigree", None)
if pedigree is None:
    raise Exception("Please provide a pedigree input.")

sample = snakemake.params.get("sample", None)
if sample is None:
    raise Exception("Please provide a sample name.")

output_dir = os.path.dirname(snakemake.output.get("vcf"))

shell("""
    /opt/edico/bin/dragen -r {ref_dir} {extra} \
        {fasta} \
        --output-directory {output_dir} \
        --output-file-prefix {sample} \
        --enable-map-align false \
        --enable-map-align-output false \
        --enable-sv true \
        --sv-denovo-scoring true \
        {bams_list} \
        --pedigree-file {pedigree} \
        --lic-server {license_txt} {log}
""")

##
# Originally also:
# --RGID {sample} \
# --RGSM {sample} \
# was included in the command line, but "when using bam input, read groups may
# not be specified on the command line"
