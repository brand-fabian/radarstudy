__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

input_bam = snakemake.input.get("bam", None)
input_bai = snakemake.input.get("bai", None)

if input_bam is None or input_bai is None:
    raise Exception("Error: Missing .bam or .bai for stats computation.")

output_stats = snakemake.output.get("stats", None)
output_flagstats = snakemake.output.get("flagstat", None)
output_idxstats = snakemake.output.get("idxstats", None)

if output_stats is None:
    raise Exception("Error: Stats needs at least an stats output location.")

extra = snakemake.params.get("extra", "")

extra_stats = extra
if not "-d" in extra_stats:
    extra_stats += " -d"

shell(r"""
    samtools stats {extra_stats} {input_bam} > {output_stats} {log}
""")
if output_flagstats is not None:
    shell(r"""
        samtools flagstat {extra} {input_bam} > {output_flagstats} {log}
    """)
if output_idxstats is not None:
    shell(r"""
        samtools idxstats {extra} {input_bam} > {output_idxstats} {log}
    """)