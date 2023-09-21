import pytest
import snakemake
import os
import shutil
import re

output_files = {
    "md.bam": "e9e2db8558cb314d430378ade747e3b6",
    "md.bai": "98097aa2252c9dc864d1faf3f81565f4",
    "md.bam.md5": "76d0c89f9a1bb12d1446ce01f52abf58",
    "md.metrics.txt": "354391abedc8bae27d2e848a099ce60a",
}

def cleanup():
    for f in output_files.keys():
        os.remove(f)
    shutil.remove(".snakemake")

def test_mark_duplicates(path, cleaned_md5, capsys):
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
