from setuptools import setup
import codecs
import os.path


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


setup(
    name="evokerlite",
    version=get_version("evokerlite/evokerlite.py"),
    author="Daniel Rice",
    author_email="dlrice@ebi.ac.uk",
    description=(
        "Produce cluster plots for SNP array data (including UK Biobank release 2)."
    ),
    license="MIT",
    packages=["evokerlite"],
    install_requires=["numpy >= 1.0.0", "matplotlib >= 1.0.0", "scipy >= 0.16.0",],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
    ],
    entry_points={"console_scripts": ["evoker-lite=evokerlite.evokerlite:cli"]},
)
