from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Compiler import Options

Options.annotate = True

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("mcsdr", ["mcsdr.pyx"], language="c++")
    ]
)
