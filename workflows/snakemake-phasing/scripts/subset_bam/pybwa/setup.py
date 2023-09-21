from setuptools import Extension, setup
from Cython.Build import cythonize
from glob import glob
import pysam

sources = set(glob("bwa/*.c")) - set(["bwa/main.c"])
extensions = [
    Extension("bwa", ["bwa.pyx", *sources], libraries=["z"],
              include_dirs=["bwa/", *pysam.get_include()],
              extra_link_args=pysam.get_libraries(),
              define_macros=[("PYSAM_STREAM_H", "1")],
              ),
]

setup(
    ext_modules = cythonize(extensions, language_level=3, annotate=True)
)