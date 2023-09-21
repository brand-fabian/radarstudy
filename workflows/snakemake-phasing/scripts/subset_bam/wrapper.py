###############################################################################
# Snakemake wrapper for subset_bam script.                                    #
#                                                                             #
# This wrapper is meant to re-format a set of bam files to concur with        #
# whatsHaps data input requirements. In particular, this means merging the    #
# bam files of all members of a pedigree and setting @RG tags accordingly.    #
#                                                                             #
# This wrapper takes care of the snakemake plumbing for hooking up the script #
# in this directory.                                                          #
###############################################################################
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import os

shell.executable('/bin/bash')
current_dir = os.path.dirname(os.path.realpath(__file__))

bam_input = ""
for sample, in_bam in snakemake.input.items():
    bam_input += " {}:{}".format(sample, os.path.abspath(in_bam))

ped = snakemake.params.get("pedigree", None)
ped = "" if ped is None else "--fam " + os.path.abspath(ped)

sites = snakemake.params.get("sites", None)
if sites is None:
    raise Exception("no sites table provided")
sites = os.path.abspath(sites)

family_id = snakemake.params.get("family_id", None)
family_id = "" if family_id is None else "--family-id " + family_id

sample_id = snakemake.params.get("sample", None)
sample_id = "" if sample_id is None else "--sample-id " + sample_id 

bam_output = snakemake.output.get("bam", None)
if bam_output is None:
    raise Exception("no bam output file path provided")
bam_output = os.path.abspath(bam_output)

threads = snakemake.threads
reference = snakemake.params.get("reference", "grch37")
flanking_region = snakemake.params.get("flanking_region", 5000)
grch37 = snakemake.params.get("grch37", None)
grch38 = snakemake.params.get("grch38", None)
fetch_unmapped_pairs = snakemake.params.get("fetch_unmapped_pairs", False)
ignore_unmatched_reads = snakemake.params.get("ignore_unmatched_reads", True)

if grch37 is None and grch38 is None:
    raise Exception("please supply either hg19 or hg38 reference genome")

script = os.path.join(current_dir, 'subset_bam.py')
log = "2>{}".format(os.path.abspath(snakemake.log[0])) if snakemake.log is not None and snakemake.log != "" else ""

grch37_str = "--grch37 {}".format(os.path.abspath(grch37)) if grch37 is not None else ""
grch38_str = "--grch38 {}".format(os.path.abspath(grch38)) if grch38 is not None else ""
fetch_unmapped_pairs_str = "--fetch-unmapped-pairs" if fetch_unmapped_pairs else ""
ignore_unmatched_reads_str = "" if ignore_unmatched_reads else "--no-ignore-unmatched-reads"
flanking_region_str = "--flanking-region {}".format(flanking_region)
threads_str = "--threads {}".format(threads) if threads is not None and threads > 1 else ""

os.makedirs(os.path.dirname(bam_output), exist_ok=True)

shell(dedent("""
    set -x
    python3 {script} --verbose debug \
        -S {sites} --reference {reference} {grch37_str} {grch38_str} {flanking_region_str} \
        {fetch_unmapped_pairs_str} {ignore_unmatched_reads_str} {threads_str} \
        {family_id} {sample_id} \
        {ped} \
        {bam_input} {log} | samtools sort -O bam -o {bam_output}
    samtools index {bam_output}
"""))