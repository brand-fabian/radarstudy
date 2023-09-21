__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent
from collections.abc import Iterable

input_bcf = snakemake.input.get("bcf", None)
if input_bcf is None:
    raise Exception("Please provide at least one input bcf file.")
if isinstance(input_bcf, str):
    input_bcf = [ input_bcf ]
if not isinstance(input_bcf, Iterable):
    raise Exception("Please provide a list of bcf files as input.")
input_bcf = list(map(os.path.abspath, input_bcf))
if any(not os.path.isfile(f) for f in input_bcf):
    raise Exception("Some bcf files appear to be missing.")
bcf_str = ' '.join(input_bcf)

output = snakemake.output

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
module_name = snakemake.params.get("module_name", "delly")

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
    delly merge {extra} -o {output} {bcf_str}
"""))