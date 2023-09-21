import pytest
import snakemake
import os
import re
import shutil

def cleanup():
    os.remove("metrics.txt")
    shutil.rmtree(".snakemake")

def test_collect_hs_metrics(path, cleaned_md5, capsys):
    result = snakemake.snakemake(
        "Snakefile",
        use_conda=True
    )
    assert result
    assert os.path.isfile("metrics.txt")
    assert cleaned_md5("metrics.txt") == "d9e5e5eb732403cd75674f2b65938b40"
    output = capsys.readouterr()
    assert re.search(r"Finished job 0", output.err)
    cleanup()