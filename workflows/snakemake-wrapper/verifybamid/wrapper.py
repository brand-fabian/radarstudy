##
# Snakemake wrapper script for gatk SelectVariants
# Caution: This wrapper is using GATK v4.1
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
import uuid
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
sample_name = snakemake.params.get("sample_name", None)

known_suffix = {
    "self_rg": ".selfRG",
    "self_sm": ".selfSM",
    "best_rg": ".bestRG",
    "best_sm": ".bestSM",
    "depth_rg": ".depthRG",
    "depth_sm": ".depthSM"
}

prefix = snakemake.params.get("prefix", str(uuid.uuid4()))

base = ".".join(snakemake.input.bam.split(".")[:-1])

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    output_prefix = os.path.join(tempdir, prefix)
    shell(r"""
        echo "Target prefix: {output_prefix}"
        verifyBamID {extra} \
            --vcf {snakemake.input.vcf} \
            --bam {snakemake.input.bam} \
            --out {output_prefix} {log}
    """)

    ##
    # Move desired output files to target location
    for key, value in known_suffix.items():
        if key in snakemake.output.keys():
            source_path = output_prefix + value
            target_path = snakemake.output.get(key)
            if target_path != source_path:
                if sample_name is not None:
                    ##
                    # Replace first column of file with sample_name (to avoid
                    # read group overlaps in subsequent analysis)
                    shell(r"""
                        perl -p -i -e \"s/^[^\\#]\(.*?\)\\t/{sample_name}/g\" {source_path}
                    """)
                shell("mv {source_path} {target_path}")

