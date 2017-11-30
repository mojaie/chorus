from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Compiler import Options

Options.annotate = True

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension(
            "mcsdr", sources=["mcsdr.pyx"], language="c++",
            extra_compile_args=["-std=c++11", "-stdlib=libc++"],
            extra_link_args=["-std=c++11", "-stdlib=libc++"]
        )
    ]
)
