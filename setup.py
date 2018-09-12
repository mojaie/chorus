
from codecs import open
import json
from os import path
from setuptools import setup, find_packages, Extension
import sys


here = path.abspath(path.dirname(__file__))

with open(path.join(here, "config.json"), "rt") as f:
    conf = json.load(f)

# Generate meta.yaml
metayaml = conf["metayaml"]
metayaml["package"]["version"] = conf["setuppy"]["version"]
metayaml["source"]["git_rev"] = conf["setuppy"]["version"]
metayaml["about"]["home"] = conf["setuppy"]["url"]
metayaml["about"]["license"] = conf["setuppy"]["license"]

with open(path.join(here, "./conda-recipe/meta.yaml"), "wt") as f:
    json.dump(metayaml, f)

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

# setup
setup_dict = conf["setuppy"]
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    setup_dict["long_description"] = f.read()

setup_dict["packages"] = find_packages(exclude=["chorus.test*"])
setup_dict["package_data"] = {
    "": ["*.yaml", "resources/test/*.mol", "resources/DrugBank/*.mol"]
}
setup_dict["ext_modules"] = [ext]
setup_dict["cmdclass"] = {'build_ext': build_ext}

setup(**setup_dict)
