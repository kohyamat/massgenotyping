import codecs
import os.path

from setuptools import find_packages, setup


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="massgenotyping",
    version=get_version("massgenotyping/__init__.py"),
    description=(
        "Python package for microsatellite genotyping from amplicon sequencing data"
    ),
    long_description=long_description,
    url="https://github.com/kohyamat/massgenotyping",
    author="Tetsuo I. Kohyama",
    author_email="tetsuo.kohyama@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    keywords="genotyping microsatellite NGS amplicon-sequencing",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.74",
        "dataclasses;python_version=='3.6'",
        "dataclasses_json",
        "fuzzysearch>=0.6.2",
        "matplotlib>=3.0.3",
        "natsort>=5.1.0",
        "numpy>=1.16.2",
        "python-Levenshtein>=0.12.0",
        "seaborn>=0.5.0",
        "tqdm>=4.30.0"
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "mgt = massgenotyping:main",
        ]
    },
)
