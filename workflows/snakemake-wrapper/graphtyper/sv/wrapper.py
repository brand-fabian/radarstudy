__author__ = "Fabian Brand"
__email__ = "brand@imbie.uni-bonn.de"

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
import os


shell.executable("/bin/bash")


bams = snakemake.input.get("bams", [])
if len(bams) == 0:
    raise Exception("Please provide at least one input bam")

ref = snakemake.params.get("reference", None)
if ref is None:
    raise Exception("No reference genome found")

vcf_in = snakemake.input.get("vcf", None)
if vcf_in is None:
    raise Exception("Missing input vcf path.")

regions = snakemake.input.get("regions", None)
regions_str = "" if regions is None else "--region_file {}".format(regions)

vcf_out = snakemake.output.get("vcf", None)
if vcf_out is None:
    raise Exception("Missing VCF output path")

extra = snakemake.params.get("extra", "")
chromosomes = snakemake.params.get("chromosomes", [*map(str, range(1, 23)), "X", "Y"])
module_name = snakemake.params.get("module_name", "bio/graphtyper")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads or 1

if len(bams) < threads:
    print(
        "WARNING: Graphtyper only supports threads up to the number of bam/cram files in the input."
    )


with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    bcftools_cmd = "-Ob" if "bcf" in vcf_out else "-Oz"
    bcftools_cmd += " -o {vcf_out}".format(vcf_out=vcf_out)
    bams = "\n".join(bams)
    shell(
        dedent(
            """
        set -ex
        
        # Load graphtyper module if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module use /software/easybuild/modules/all
            module load {module_name}
        fi

        export TMPDIR="{tempdir}"
        export TMPDIR=$(mktemp -d)

        {{
            # Create sams list
            cat <<EOS >{tempdir}/sam.list
            {bams}
            EOS

            # Create output dir
            mkdir -p {tempdir}/output

            graphtyper genotype_sv {ref} {vcf_in} {extra} --sams {tempdir}/sam.list \
                {regions_str} \
                --threads {threads} \
                --output {tempdir}/output

            # Concat the output files of graphtyper into the results
            echo {chromosomes} | tr ' ' '\\n' | while read chrom; do
                if [[ ! -d {tempdir}/output/${{chrom}} ]]; then
                    continue
                fi
                find {tempdir}/output/${{chrom}} -name "*.vcf.gz" | sort
            done >{tempdir}/output.vcf.list
            bcftools concat --naive --file-list {tempdir}/output.vcf.list | bcftools sort {bcftools_cmd}
            tabix {vcf_out}
        }} {log}
            """
        )
    )
