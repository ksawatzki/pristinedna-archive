# PristineDNA

An R package to process environmental DNA (eDNA) data and format it for submission to OBIS and GBIF.

## Installation

You can install this package from GitHub:

```r
# install.packages("devtools")  # If not installed
devtools::install_github("yourusername/PristineDNA")
```

## Usage

```r
library(PristineDNA)

PristineDNA(metadata = "metadata.xlsx", eflow = "eDNA_data.tsv", seq = "sequences.fasta")
```

