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
  metadata = "metadata.xlsx",      # Lab spreadsheet with metadata and library prep records
  eflow    = "eDNA_data.tsv",      # TSV from your eDNA pipeline
  seq      = "sequences.fasta",    # optional FASTA with sample IDs in headers
  qc       = TRUE,                 # optional QC with screen and printed output
)

# NCBI SRA spreadsheets
pristinesra(
  metadata       = "metadata.xlsx",         # Lab spreadsheet with metadata and library prep records
  outdir         = "sra_out",               # where TSVs are written
)

pristinedna(metadata = "metadata.xlsx", eflow = "eDNA_output.tsv", qc = TRUE)
pristinesra(metadata = "metadata.xlsx")
```
## Output

```pgsql
dwca_out/
  ├─ occurrence.txt
  ├─ dna.txt
  ├─ eml.xml
  ├─ meta.xml
  └─ project_dwca.zip
  ├─ qc_summary_occurrence.txt
  ├─ qc_summary_dnaderiveddata.txt


sra_out/
  ├─ BioSample.tsv
  ├─ SRAmetadata.tsv
```

## Input

1) metadata.xlsx (pristinedna, pristinesra)
  - Project/Events sheet (project-level + event-level)
    - Minimum required columns: ps_sample_id (join key), date, lat, long,
    collection_time, country, env_broad_scale, env_medium.
  - Samples sheet (sample-level)
    - Minimum required columns: ps_sample_id (join key)
    - Will prompt for additional required data if blank
  - Specific output may require additional columns

2) eDNA_data.tsv (pristinedna)
  - Minimum required columns: ps_sample_id column names (to link to metadata),
  minimum of 1 taxonomic idenitifying unit (no taxa OTUs will be dropped), pid
  (percent identity to assigned taxa)
  - OTU (sequence unit ID)
  - genus and species
  - read counts per sample-OTU
