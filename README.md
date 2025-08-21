# PristineDNA

An R package to standardize environmental DNA (eDNA) outputs and package them for submission to **OBIS/GBIF** (Darwin Core Archive, DwC-A) and **NCBI SRA** from the same project metadata.
The package follows Pristine Seas metadata conventions and MIxS/MIMARKS-style fields.

- `pristinedna()` — builds a DwC-A (zip) for OBIS/GBIF  
- `pristinesra()` — builds submission tables for NCBI SRA

---

## Key features

- **Single entry points**: one-call workflows to quickly format eDNA submissions to multiple repositories.
- **Validation**: checks presence of required fields and formats before writing.
- **Darwin Core Archive export**: writes a zipped archive plus component files (`occurrence.txt`, `dna.txt`, `eml.xml`, `meta.xml`).
- **SRA builder**: generates biosample TSVs for direct upload to NCBI repositories.

---

## Installation

```r
install.packages("devtools")
devtools::install_github("ksawatzki/PristineDNA")
# or: remotes::install_github("ksawatzki/PristineDNA")
```

## Quick start

```r
library(PristineDNA)

# OBIS/GBIF (DwC-A)
pristinedna(
  metadata = "metadata.xlsx",      # Excel workbook with 2 sheets (see below)
  eflow    = "eDNA_data.tsv",      # TSV from your eDNA pipeline (long format)
  seq      = "sequences.fasta",    # optional FASTA with ZOTU/ASV IDs in headers
  outdir   = "dwca_out",           # optional; default = working directory
  archive  = "project_dwca.zip"    # optional; default = auto-named
)

# NCBI SRA spreadsheets
pristinesra(
  metadata       = "metadata.xlsx",         # Lab spreadsheet with metadata and library prep records
  outdir         = "sra_out",               # where TSVs/XML are written
  platform       = "ILLUMINA",              # e.g., ILLUMINA, OXFORD_NANOPORE, PACBIO_SMRT
  library_prep   = list(strategy="AMPLICON", selection="PCR", layout="PAIRED"),
  validate       = TRUE
)

pristinedna(metadata = "metadata.xlsx", eflow = "eDNA_output.tsv", seq = "seqs.fasta")
pristinesra(metadata = "metadata.xlsx", fastq_manifest = man, outdir = tempdir())
```
## Output

```pgsql
dwca_out/
  ├─ occurrence.txt
  ├─ dna.txt
  ├─ eml.xml
  ├─ meta.xml
  └─ project_dwca.zip

sra_out/
  ├─ BioSample.tsv
  ├─ SRAmetadata.tsv
  └─ SRA_Submission.xml
```

## Input

1) metadata.xlsx (pristinedna, pristinesra)
  - Project/Events sheet (project-level + event-level)
    - Recommended: ps_sample_id (join key), expedition, depth_m,
    decimalLatitude, decimalLongitude, country, env_broad_scale, env_medium.
  - Samples sheet (sample-level)
    - Recommended: ps_sample_id (join key), lib_layout, seq_meth, marker.

2) eDNA_data.tsv (pristinedna)
  - ps_sample_id column names (to link to metadata)
  - asv / zotu (sequence unit ID)
  - genus and species (assigned taxon)
  - readCount (abundance proxy)
