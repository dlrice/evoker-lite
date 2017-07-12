from setuptools import setup

setup(
    name = "evoker_lite",
    version = "0.0.2",
    author = "Daniel Rice",
    author_email = "dr9@sanger.ac.uk",
    description = ("Produce cluster plots for SNP array data (including UK Biobank release 2)."),
    license = "MIT",
    packages=["evokerlite"],
    install_requires=['numpy >= 1.0.0',
                      'matplotlib >= 1.0.0',
                      'scipy >= 0.16.0',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
    ],
)