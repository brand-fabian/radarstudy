__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell

shell.executable('/bin/bash')

tool = "bedtools intersect"

has_log = snakemake.log is None or snakemake == ""
has_params = snakemake.params is None or snakemake.params == "" or snakemake.params.intersect is None or snakemake.params.intersect == ""

if has_log and has_params:
    shell(tool + " {snakemake.params.intersect} -a {snakemake.input.a} -b {snakemake.input.b} > {snakemake.output[0]} 2> {snakemake.log}")
elif has_log:
    shell(tool + " {snakemake.input} > {snakemake.output[0]} 2> {snakemake.log}")
elif has_params:
    shell(tool + " {snakemake.params.intersect} -a {snakemake.input.a} -b {snakemake.input.b} > {snakemake.output[0]}")
else:
    shell(tool + " -a {snakemake.input.a} -b {snakemake.input.b} > {snakemake.output[0]}")