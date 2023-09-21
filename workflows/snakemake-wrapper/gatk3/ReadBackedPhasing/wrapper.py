##
# Snakemake wrapper script for gatk haplotype caller
# GATK: This wrapper is using GATK v3.8, install it into the environment
#       prior to usage
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")

memory = snakemake.resources.get("mem_mb", None)
if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(r"""
    gatk3 {jvm_args} -T ReadBackedPhasing \
        -R {snakemake.input.ref} \
        -I {snakemake.input.bam} \
        -V {snakemake.input.vcf} \
        -L {snakemake.input.vcf}
        -o {snakemake.output.vcf} \
        {extra} {log}
""")