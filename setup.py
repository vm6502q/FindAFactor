import os
import setuptools
from distutils.core import setup, Extension

cpp_args = ['-std=c++17', '-lpthread']

README_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.md')
with open(README_PATH) as readme_file:
    README = readme_file.read()

ext_modules = [
    Extension(
        '_find_a_factor',
        ["FindAFactor/_find_a_factor.cpp", "FindAFactor/dispatchqueue.cpp"],
        include_dirs=['FindAFactor/include', 'pybind11/include'],
        language='c++',
        extra_compile_args = cpp_args,
    ),
]

setup(
    name='FindAFactor',
    version='1.0.3',
    author='Dan Strano',
    author_email='dan@unitary.fund',
    description='Find any nontrivial factor of a number',
    long_description=README,
    long_description_content_type='text/markdown',
    url="https://github.com/vm6502q/FindAFactor",
    license="MIT",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=["pybind11"],
    ext_modules=ext_modules,
    packages=['FindAFactor'],
)
