##
# Snakemake wrapper script for ngs-chew fingerprint
# 
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
import sys
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
reference = snakemake.input.get("ref", None)
bam = snakemake.input.get("bam", None)
output_fingerprint = snakemake.output.get("fingerprint", snakemake.output[0])
output_vcf = snakemake.output.get("vcf", None)

if reference is None:
    raise Exception("Error: No reference build provided.")

if bam is None:
    raise Exception("Error: No bam provided.")

sample_name = bam.split("/")[-1].split(".")[0]

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    print("Using temp dir {}".format(tempdir), file=sys.stderr)
    output_prefix = "{}/{}".format(tempdir, sample_name)
    shell(r"""
        ngs-chew fingerprint \
            --reference {reference} \
            --output-fingerprint {output_prefix} \
            --input-bam {bam} \
            --write-vcf
    """)
    shell(r"""
        mv {output_prefix}.npz {output_fingerprint}
    """)

    if output_vcf is not None:
        shell(r"""
            mv {output_prefix}.vcf.gz {output_vcf}
            tabix {output_vcf}
        """)