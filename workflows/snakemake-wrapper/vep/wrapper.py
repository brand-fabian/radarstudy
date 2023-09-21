__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

input_vcf = snakemake.input.get('vcf', None)
if input_vcf is None:
    raise Exception("Please provide an input vcf.")

output_vcf = snakemake.output.get('vcf', None)
if output_vcf is None:
    raise Exception("Please provide an output vcf path")

text_output = 'vcf' not in output_vcf
if text_output:
    output_options = ''
else:
    if output_vcf.endswith('gz'):
        output_options = '--vcf --compress_output bgzip'
    else:
        output_options = '--vcf'

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
##
# Manage run in temporary directory
#
with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    cache_dir = snakemake.params.get('vep_cache', None)
    reference_build = snakemake.params.get('reference_build', 'GRCh37')
    species = snakemake.params.get('species', 'homo_sapiens')
    if cache_dir is None or not os.path.exists(cache_dir):
        # Download cache
        shell("""
            vep_install -a cf -c {tempdir} -s {species} -y {reference_build} -t {log}
        """)
        log = log_append
        cache_dir = tempdir
    
    extra = snakemake.params.get('extra', '')
    shell("""
        vep --force_overwrite {extra} -a {reference_build} --cache --dir {cache_dir} -i {input_vcf} {output_options} -o {output_vcf} {log}
    """)
    log = log_append
    if not text_output:
        shell("""
            tabix {output_vcf} {log}
        """)