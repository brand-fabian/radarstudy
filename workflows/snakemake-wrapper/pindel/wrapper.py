# Pindel Snakemake Wrapper 
__author__ = 'Sugi Sivalingam'
__email__ = 's.sivalingam@uni-bonn.de'

import os
import uuid
from snakemake.shell import shell
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

# pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
# D = Deletion, SI = Short Insertion, INV = Inversion, TD = Tandem Duplication, LI = Large Insertion, BP = Unassigned Breakpoints, RP= Read Pair

#known_suffix = {
#    "D": "_D",
#    "BP": "_BP",
#    "INV": "_INV",
#    "TD": "_TD",
#    "LI": "_LI",
#    "SI": "_SI"
#    "RP": "_RP"
#}

prefix = snakemake.params.get("prefix", str(uuid.uuid4()))
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads if hasattr(snakemake, 'threads') else 1

input_config = snakemake.input.get("config", None)
if input_config is None:
    raise Exception("No input config provided. Please provide at least [bam]")

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    output_prefix = os.path.join(tempdir, prefix)
    shell(r"""
        pindel -T {snakemake.threads} {snakemake.params.extra} -f {snakemake.input.ref} -i {input_config} -o {output_prefix} {log}
        pindel2vcf -P {output_prefix}  -r {snakemake.input.ref}  -R "human_g1k_v37" -d "human_g1k_v37" -v {snakemake.output.vcf}
        bcftools filter 
    """)

