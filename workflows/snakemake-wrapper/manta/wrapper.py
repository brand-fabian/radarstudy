__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from uuid import uuid4
from tempfile import mkdtemp, gettempdir
from collections.abc import Iterable
import textwrap
import os
import shutil


shell.executable('/bin/bash')

output_files = {
    "diploid": "diploidSV.vcf.gz",
    "candidate": "candidateSV.vcf.gz",
    "candidateSmallIndels": "candidateSmallIndels.vcf.gz",
}

def process_input(data, name):
    if data is None:
        raise Exception("No input provided. Please provide at least one {}.".format(name))
    if isinstance(data, str):
        data = [ data ]
    if not isinstance(data, Iterable):
        raise Exception("{} must be an iterable of bam file paths.".format(name))
    return data

input_bam = process_input(snakemake.input.get("bam", None), "bam")

call_denovo = True if snakemake.params.get("call_denovo", False) else False
pedigree = process_input(snakemake.params.get("pedigree", None), "pedigree")
if call_denovo and len(pedigree) != 3:
    raise Exception("Please provide a pedigree array of form [<child_id> <father_id> <mother_id>]; Got {}".format(pedigree))

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide the reference fasta as parameter.")
reference = os.path.abspath(reference)

prefix = snakemake.params.get("prefix", str(uuid4()))
extra_config = snakemake.params.get("extra_config", "")
extra_caller = snakemake.params.get("extra_caller", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_caller = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
memory = snakemake.resources.get("mem_mb", 1024 * 1024 * 4) // 1024

workdir = mkdtemp(dir=os.getenv("SCRATCH_DIR", gettempdir()), prefix=prefix)
try:
    bam_str = '--bam ' + ' --bam '.join(input_bam)
    shell(textwrap.dedent("""
        set -ex
        configManta.py {extra_config} {bam_str} \
            --referenceFasta {reference} \
            --runDir {workdir} {log}

        python {workdir}/runWorkflow.py -j {threads} -g {memory} {log_caller}
    """))

    if call_denovo:
        ped_str = ' '.join(pedigree)
        shell(textwrap.dedent("""
            set -ex
            $(dirname $(which configManta.py))/../share/manta-1.6.0-1/libexec/denovo_scoring.py \
                {workdir}/results/variants/diploidSV.vcf.gz {pedigree}
        """))

    for key, value in snakemake.output.items():
        if key in output_files:
            if not os.path.isdir(os.path.dirname(value)):
                os.makedirs(os.path.dirname(value), exist_ok=True)

            shell(textwrap.dedent(f"""
                set -ex
                mv -v {workdir}/results/variants/{output_files[key]} {value}
                tabix {value}
            """))
finally:
    shutil.rmtree(workdir)