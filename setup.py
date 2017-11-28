
from codecs import open
from os import path
from setuptools import setup, find_packages, Extension
import sys

ignore_cython = False
ext = None

try:
    from Cython.Distutils import build_ext
    ext = Extension(
        "chorus.cython.mcsdr", ["chorus/cython/mcsdr.pyx"], language="c++"
    )
except ImportError:
    if ignore_cython:
        from distutils.command.build_ext import build_ext
        ext = Extension(
            "chorus.cython.mcsdr", ["chorus/cython/mcsdr.cpp"]
        )
        print("Cython is not available. Distribute C++ source files.")
    else:
        print("Error: Cython is required.")
        sys.exit()

here = path.abspath(path.dirname(__file__))
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="Chorus",
    version="0.7.2",
    description="Simple chemical structure modeling toolkit",
    long_description=long_description,
    url="https://github.com/mojaie/chorus",
    author="Seiji Matsuoka",
    author_email="mojaie@aol.com",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6"
    ],
    keywords="molecule-modeling cheminformatics",
    packages=find_packages(exclude=["chorus.test*"]),
    package_data={
        "": ["*.yaml", "resources/test/*.mol", "resources/DrugBank/*.mol"]
    },
    python_requires=">=3.5",
    install_requires=[
        "numpy", "numexpr", "matplotlib", "networkx", "pyyaml"
    ],
    ext_modules=[ext],
    cmdclass={'build_ext': build_ext}
)
