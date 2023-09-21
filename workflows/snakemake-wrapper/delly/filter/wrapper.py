__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent

input_bcf = snakemake.input.get("bcf", None)
if input_bcf is None:
    raise Exception("Please provide an input bcf file.")
input_bcf = os.path.abspath(input_bcf)
if not os.path.isfile(input_bcf):
    raise Exception("The bcf file could not be found at {}.".format(input_bcf))

output = snakemake.output

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
module_name = snakemake.params.get("module_name", "delly")

filter_type = snakemake.params.get("filter", "germline")
if filter_type not in ["germline", "somatic"]:
    raise Exception("Delly filter type must be set. Valid Choices: [ germline, somatic ]")

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
    export OMP_NUM_THREADS={threads}
    delly filter {extra} -f {filter_type} -o {output} {input_bcf}
"""))