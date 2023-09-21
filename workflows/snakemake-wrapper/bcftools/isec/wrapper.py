__author__ = "Fabian Brand"
__email__ = "brand@imbie.meb.uni-bonn.de"

from snakemake import shell
import tempfile
import os

shell.executable("/bin/bash")

exclude = snakemake.params.get('exclude', None)
include = snakemake.params.get('include_expr', None)
isec = snakemake.params.get('isec', None)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if not ('-T' in isec or '-t' in isec) and len(snakemake.input.variants) < 2:
    raise Exception("One variants file requires --targets option for bcftools.")

with tempfile.TemporaryDirectory() as tmpdir:
    command = "bcftools isec"
    if exclude is not None:
        command += " -e \"{}\"".format(exclude)
        del snakemake.params.exclude
    if include is not None:
        command += " -i \"{}\"".format(include)
        del snakemake.params.include_expr
    if len(snakemake.input) == 1:
        command += " -o {}".format(snakemake.output.variants)
    else:
        command += " -p {}".format(tmpdir)
    if snakemake.output[0].endswith(".gz"):
        command += " -O z"
    if isec is not None:
        command += " {}".format(isec)
    for i in snakemake.input:
      command += " {}".format(i)
    if snakemake.log:
        command += " {}".format(log)
    if len(snakemake.input) > 1:
        for i, out in enumerate(snakemake.output):
            if out.endswith('.gz'):
                command += " && cp {}/{:04d}.vcf.gz {}".format(tmpdir, i, out)
            else:
                command += " && cp {}/{:04d}.vcf {}".format(tmpdir, i, out)
    shell(command)
