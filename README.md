# Evoker Lite

Evoker Lite is a python tool to generate cluster plots (PNG/PDF files) given a set of PLINK + intensity files. It is Lite in the sense that the original [Evoker](https://github.com/wtsi-medical-genomics/evoker) is a Java program which allows for interactive plotting (zooming, panning) and also selection and reassignment of calls.

## UK Biobank v2

The current release of this tool has been designed to work with v2 of the UK Biobank data. For example, to generate cluster plots for each of the rsid's named in `rsid_list.txt` and over all of the batches:

 ```
evoker-lite \
 --ukb \
 --data /ukbiobank/release/Genotypes/ \
 --fam ~/ukb1234_cal_v2_s488377.fam \
 --output ~/plots/ \
 --rsids ~/rsid_list.txt
 ```

 Once complete the directory `~/plots/`, in this example, will contain subdirectories for each rsid and within these a PNG for each batch.


**Note** The original [Evoker](https://github.com/wtsi-medical-genomics/evoker) is presently being developed to allow viewing/call reassignment of UK Biobank data.

## Command line arguments

```
--ukb               flag to indicate data is in UKBiobank format
-d, --data          directory of PLINK/intensity data
-f, --fam           location of the fam file if not in the directory specified with -d/--data
-o, --output        directory to save plots
-r RSIDS, --rsids   text file with rsids to create cluster plots from
--no-transform      [UK Biobank] flag to not use contrast/strength coordinates
--snp-posterior     [UK Biobank] flag to plot SNP Posterior
```

## Requirements

 * Python 2/3
 * numpy
 * scipy
 * matplotlib

The quickest and easiest way to install all of these is via [miniconda](https://conda.io/miniconda.html) or [anaconda](https://www.continuum.io/downloads).

## Installation

```
pip install git+git://github.com/dlrice/evoker-lite.git --user
```

