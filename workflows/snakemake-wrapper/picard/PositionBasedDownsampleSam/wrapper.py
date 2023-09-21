__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
import os

shell.executable('/bin/bash')

bam = snakemake.input.get('bam', None)
if bam is None:
    raise Exception("No input provided")
bam = os.path.abspath(bam)

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref = os.path.abspath(ref)

bam_out = snakemake.output.get("bam", None)
if bam_out is None:
    raise Exception("Please provide an output bam path")
bai_out = snakemake.output.get("bai", None)
md5_out = snakemake.output.get("md5", None)

threads = snakemake.threads
factor = snakemake.params.get("factor", 1)
if factor == 1:
    raise Exception("please set a factor 0 < factor < 1 for downsampling")

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

create_index_str = "--CREATE_INDEX" if bai_out is not None else ""
create_md5_str = "--CREATE-MD5-FILE" if md5_out is not None else ""

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex
        export TMPDIR={tempdir}

        picard {jvm_args} PositionBasedDownsampleSam {extra} \
            --FRACTION {fraction} \
            --REFERENCE {reference} \
            --INPUT {bam} \
            --OUTPUT {bam_out} {log}
    """))
