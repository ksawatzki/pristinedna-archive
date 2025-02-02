# PristineDNA

An R package to process environmental DNA (eDNA) data and format it for submission to OBIS and GBIF. This package is intended to conform to Pristine Seas standard metadata (and output files and is not directly applicable to all input formats.

## Installation

You can install this package from GitHub:
```r
install.packages("devtools")

devtools::install_github("ksawatzki/PristineDNA")
```

## Usage
```r
library(PristineDNA)

PristineDNA(metadata = "metadata.xlsx", eflow = "eDNA_data.tsv", seq = "sequences.fasta")
```

## Input
metadata are Excel spreadsheets with two sheets ...

eflow are standard tsv file output of eFlow pipeline

seq are currently optional, expectation is a FASTA file with ZOTU sample names

## Output
Program currently outputs a zipped DwC-a format file that is directly submittable to OBIS and GBIF

It also writes the four files in this zip for easy viewing. These are the occurrence file, the DNA extension, eml and meta files.

Future output will be the associated FASTA file with modified titles to include metadata. This output will be made for direct submission to GenBank or other repositories.
