import pytest
import os
import hashlib
import tempfile

@pytest.fixture
def path(request):
    """Change the current working directory into the dir containing the test
       file, s.t. snakemake runs can start from there."""
    cur_dir = os.getcwd()
    os.chdir(os.path.dirname(request.fspath))
    yield (request.fspath, request.function)
    os.chdir(cur_dir)

@pytest.fixture
def md5():
    """Compute the md5 hash of the file at the given file pointer"""
    def _md5(fp):
        hash_md5 = hashlib.md5()
        for chunk in iter(lambda: fp.read(4096), b""):
            hash_md5.update(chunk)
        return hash_md5.hexdigest()
    return _md5


@pytest.fixture
def cleaned_md5(md5):
    def _md5(f_name, binary=False):
        assert os.path.isfile(f_name)
        if binary:
            return md5(open(f_name, 'rb'))
        else:
            with tempfile.TemporaryFile() as out_f:
                with open(f_name, 'r') as in_f:
                    for line in in_f:
                        if not line.startswith("#"):
                            out_f.write(str.encode(line))
                out_f.seek(0)
                return md5(out_f)
    return _md5