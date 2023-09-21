__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent
from collections.abc import Iterable
from tempfile import TemporaryDirectory, gettempdir

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("No bam file specified.")
input_bam = os.path.abspath(input_bam)

sample = snakemake.params.get("sample", input_bam.split(".")[0])
reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide a path to a reference genome.")
reference = os.path.abspath(reference)

output_files = {
    "vcf": "{}-smoove.genotyped.vcf.gz".format(sample)
}

extra = snakemake.params.get("extra", "")

if "--genotype" not in extra:
    extra += " --genotype"
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
module_name = snakemake.params.get("module_name", "smoove")

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
            if [ -d "/ssd" ]; then
                export SINGULARITY_BINDPATH="/ceph01,/gpfs01,/ssd,/tmp"
            else
                export SINGULARITY_BINDPATH="/ceph01,/gpfs01,/tmp"
            fi
        fi
        if [ -n "$SCRATCH_DIR" ]; then
            export TMPDIR="$SCRATCH_DIR"
            export TMPDIR=$(mktemp -d)
        else
            export TMPDIR=$(mktemp -d)
            trap "rm -rf $TMPDIR" EXIT
        fi

        smoove call --outdir {tempdir} \
            {extra} --name {sample} \
            --fasta {reference} \
            -p {threads} {input_bam} {log}
    """))

    for key, value in output_files.items():
        out = snakemake.output.get(key, None)
        if out is not None:
            if value.endswith('.gz'):
                file_type = '.'.join(value.split('.')[-2:])
            else:
                file_type = value.split('.')[-1]
            
            if file_type == 'vcf':
                convert_out = f"{tempdir}/rhdr.vcf"
                convert = f"-O v -o {convert_out}"
            elif file_type == 'bcf':
                convert_out = f"{tempdir}/rhdr.bcf"
                convert = f"-O b -o {convert_out}"
            elif file_type == 'vcf.gz':
                convert_out = f"{tempdir}/rhdr.vcf.gz"
                convert = f"-O z -o {convert_out}"
            else:
                convert_out = f"{tempdir}/{value}"
                convert = None

            if convert is not None:
                # Ensure sample names are set correctly...
                shell(dedent(f"""
                    set -ex
                    echo {sample} > {tempdir}/samples.txt
                    bcftools reheader -s {tempdir}/samples.txt {tempdir}/{value} | bcftools convert {convert}
                """))
            shell(dedent("""
                set -ex
                if [ ! -d $(dirname "{out}") ]; then
                    mkdir -p $(dirname "{out}")
                fi
                mv -v {convert_out} {out} {log_append}
            """))

            if any(out.endswith(e) for e in [ 'vcf.gz', 'bcf', 'vcf' ]):
                shell(dedent("""
                    tabix {out}
                """))