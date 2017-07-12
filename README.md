# Evoker Lite

Evoker Lite is a python tool to generate cluster plots (PNG/PDF files) given a set of PLINK + intensity files. It is Lite as the original Evoker is more fully featured Java program which allows for interactive plotting (zooming, panning) and also selection and reassignment of calls.

## UK Biobank v2

The current release of this tool has been designed to work with v2 of the UK Biobank data. For example, to generate cluster plots for each of the rsid's named in `rsid_list.txt` and over all of the batches:

 ```
 python -m evoker_lite.py  \
 --ukb \
 --data ~/ukb/data/ \
 --output ~/plots/ \
 --rsids rsid_list.txt
 ```

 Once complete the directory `~/plots/`, in this example, will contain subdirectories for each rsid and within these a PNG for each batch.


## Command line arguments

```
--ukb               flag to indicate data is in UKBiobank format
-d, --data          directory of PLINK/intensity data
-o, --output        directory to save plots
-s RSIDS, --rsids   text file with rsids to create cluster plots from
--no-transform      flag to not plot UKBiobank data in contrast/strength coordinates
--no-snp-posterior  flag to not plot UKBiobank SNP Posterior
```

## Requirements

 * Python 2/3
 * numpy
 * scipy
 * matplotlib

The quickest and easiest way to install all of these is via [miniconda](https://conda.io/miniconda.html) or [anaconda](https://www.continuum.io/downloads).

## Installation

If you already have the above dependencies met, clone this repository and you can use the code immediately:

```
git clone https://github.com/dlrice/evoker-lite
cd evoker-lite
```

If you require pip to install the dependencies run the additional step:

```
python setup.py install --user
```
