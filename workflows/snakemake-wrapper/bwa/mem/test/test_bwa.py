import pytest
import snakemake
import os
import shutil
import re

output_files = {
    "output.bam": "661319295d665b806741c30a63fd830f",
    "output.bai": "58e80830899cd042878f1feac9aecd97",
    "output.bam.md5": "60d6bcbdf064ce7826aceacf2aa6f040",
}

def cleanup():
    os.remove("output.bam")
    os.remove("output.bai")
    os.remove("output.bam.md5")
    shutil.rmtree(".snakemake")

def test_bwa_mem(path, cleaned_md5, capsys):
    result = snakemake.snakemake(
        "Snakefile",
        use_conda=True,
    )
    assert result
    for f_name, md5 in output_files.items():
        assert os.path.isfile(f_name)
        assert cleaned_md5(f_name, binary=True) == md5
    output = capsys.readouterr()
    assert re.search(r"Finished job 0", output.err)
    cleanup()