__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
import os

shell.executable('/bin/bash')


class VQSRResource:
    def __init__(self, **kwargs):
        self.name = str(kwargs['name'])
        self.known = 'true' if bool(kwargs['known']) else 'false'
        self.training = 'true' if bool(kwargs['training']) else 'false'
        self.truth = 'true' if bool(kwargs['truth']) else 'false'
        self.prior = f'{kwargs["prior"]:.1f}'
        self.file = os.path.abspath(str(kwargs['file']))
    
    def __str__(self):
        return f'--resource {self.name},known={self.known},training={self.training},truth={self.truth},prior={self.prior}:{self.file}'


def get_value_or_error(obj, key, err_str):
    val = obj.get(key, None)
    if val is None:
        raise Exception(err_str)
    return val

vcf = get_value_or_error(snakemake.input, 'vcf', 'No vcf input provided')
input_str = f"--in-vcf {vcf}"

vcf_out = get_value_or_error(snakemake.output, 'vcf', 'No vcf output provided')
recal_out = get_value_or_error(snakemake.output, 'recal', 'No recalibration table output provided')
tranches_out = get_value_or_error(snakemake.output, 'tranches', 'No tranches output provided')

compress = vcf_out.endswith('.gz')
uncompressed_vcf = '.'.join(vcf_out.split('.')[:-1])
if compress:
    output_str = f"--out-vcf {uncompressed_vcf}"
else:
    output_str = f"--out-vcf {vcf_out}"
output_str += f" --out-recal {recal_out} --out-tranches {tranches_out}"

resources = get_value_or_error(snakemake.params, 'resources', 'Please provide atleast one VQSR resource')
if not isinstance(resources, list) or len(resources) == 0:
    raise Exception('Please provide a list of resources (length >= 1)')
res_objs = list(map(lambda d: VQSRResource(**d), resources))
res_str = ' '.join(map(str, res_objs))

annotations = get_value_or_error(snakemake.params, 'annotations', 'Please provide atleast one annotation.')
if not isinstance(annotations, list) or len(annotations) == 0:
    raise Exception('Please provide a list of annotations (length >= 1)')
anno_str = "--annotation " + " --annotation ".join(annotations)

threads = snakemake.threads
extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(dedent("""
    # Setup scratch dir (local ssd/generic temp)
    SSD=${{SCRATCH_DIR:-}}
    if [ -n "$SSD" ]; then
        export TMPDIR=$SCRATCH_DIR
    fi
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    # Load parabricks if a module system is present
    if `type module 2>&1 | grep -q function`; then
        module load {module_name}
    fi

    pbrun vqsr {extra} \\
        {res_str} \\
        {input_str} \\
        {output_str} \\
        {anno_str} \\
        --tmp-dir $TMPDIR --x3 {log}
"""))

if compress:
    shell(dedent("""
        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbgzip -n {threads} {uncompressed_vcf}
        pbrun indexgvcf --input {vcf_out}
        rm -rf {uncompressed_vcf}
    """))
