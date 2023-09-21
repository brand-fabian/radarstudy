__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import uuid
from tempfile import TemporaryDirectory, gettempdir
import os
import re
from textwrap import dedent

shell.executable('/bin/bash')

known_suffix = {
    "dist": "mosdepth.global.dist.txt",
    "summary": "mosdepth.summary.txt",
    "per_base": "per-base.bed.gz",
    "regions_dist": "mosdepth.region.dist.txt",
    "regions": "regions.bed.gz",
    "quantized": "quantized.bed.gz",
    "thresholds": "thresholds.bed.gz"
}

commands = {
    "dist": lambda x: True,
    "summary": lambda x: True,
    "per_base": lambda x: re.search(r"(^|\s+)[^a-zA-Z\-]?-n", x) is None and re.search(r"(^|\s+)[^a-zA-Z\-]?--no-per-base", x) is None,
    "regions": lambda x: re.search(r"(^|\s+)[^a-zA-Z\-]?--by", x) is not None or re.search(r"(^|\s+)[^a-zA-Z\-]?-b", x) is not None,
    "regions_dist": lambda x: re.search(r"(^|\s+)[^a-zA-Z\-]?--by", x) is not None or re.search(r"(^|\s+)[^a-zA-Z\-]?-b", x) is not None,
    "quantized": lambda x: re.search(r"(^|\s+)[^a-zA-Z\-]?--quantize", x) is not None or re.search(r"(^|\s+)[^a-zA-Z\-]?-q", x) is not None,
    "thresholds": lambda x: re.search(r"(^|\s+)[^a-zA-Z\-]?--thresholds", x) is not None or re.search(r"(^|\s+)[^a-zA-Z\-]?-T", x) is not None
}

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("No input provided. Please provide at least [bam]")

output_dist = snakemake.output.get("dist", None)
output_summary = snakemake.output.get("summary", None)
if output_dist is None or output_summary is None:
    raise Exception("""Both dist and summary outputs are created.
                    Please specify target files for each.""")

prefix = snakemake.params.get("prefix", str(uuid.uuid4()))
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    output_prefix = os.path.join(tempdir, prefix)
    shell(dedent("""
        mosdepth -t {snakemake.threads} {extra} \
            {output_prefix} \
            {input_bam} {log}
    """))
    for key, value in known_suffix.items():
        if commands[key](extra):
            source_path = "{}.{}".format(output_prefix, value)
            target_path = snakemake.output.get(key, None)

            if not os.path.isdir(os.path.dirname(target_path)):
                os.makedirs(os.path.dirname(target_path), exist_ok=True)

            if target_path is None:
                raise Exception("No output file for requested analysis ({})".format(key))
            shell(dedent("""
                mv -v {source_path} {target_path}
            """))