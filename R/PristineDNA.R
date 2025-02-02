# Check and install required packages if missing
required_packages <- c("magrittr","countrycode", "measurements", "readxl", "readr", "dplyr", 
                       "stringr", "Biostrings", "reshape2", "purrr", "XML", "utils", "zip")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages, dependencies = TRUE)
}

# Load  libraries
library(magrittr) #fixing a dplyr issue
library(countrycode)  # country code conversion
library(measurements)  # lat/long conversion
library(readxl)  # input
library(readr)   # may not be necessary anymore - KS check
library(dplyr)   # data manipulation
library(stringr) # text manipulation
library(Biostrings) # fasta file
library(reshape2) # data melting 
library(purrr) # lat/long conversion
library(XML)
library(utils)
library(zip)

# Function to load eDNA data and map to occurrence and DNA-derived data output format
pristinedna <- function(metadata, eflow, seq = NULL) {
  # Read Excel file (metadata and labwork sheets)
  metadata_df <- read_excel(metadata, sheet = "metadata")
  labwork_df <- read_excel(metadata, sheet = "labwork")
  # Read TSV file (eDNA count data)
  edna_data <- read.delim(eflow, na = c("NA", "dropped"),header=TRUE)
  
  # QC Check for pid values <90
  low_pid_rows <- edna_data %>% dplyr::filter(pid < 90)
  
  if (nrow(low_pid_rows) > 0) {
    repeat {
      cat("Percent identity values <90% have been found in the eFlow input.\n")
      cat("Do you want to:\n")
      cat("1 - Ignore and continue\n")
      cat("2 - Drop all <90%\n")
      cat("3 - Decide for each instance\n")
      cat("4 - Quit program to investigate\n")
      user_choice <- readline("Enter option (1/2/3/4): ")
      
      if (user_choice == "1") {
        break
      } else if (user_choice == "2") {
        edna_data <- edna_data %>% filter(pid >= 90)
        break
      } else if (user_choice == "4") {
        stop("See you later!")
      } else if (user_choice == "3") {
        for (i in 1:nrow(low_pid_rows)) {
          print(low_pid_rows[i, c("genus", "species", "sequence", "seq_length", "OTU", "unique_hits", "pid")])
          repeat {
            decision <- readline("1 - Keep
2 - Drop
3 - Quit program
Enter choice: ")
            if (decision == "1") {
              break
            } else if (decision == "3") {
              stop("See you later!")
              break
            } else if (decision == "2") {
              edna_data <- edna_data %>% filter(!(pid == low_pid_rows$pid[i] & genus == low_pid_rows$genus[i] & species == low_pid_rows$species[i]))
              break
            } else {
              message("Invalid entry. Please select 1 or 2.")
            }
          }
        }
        break
      } else {
        message("Invalid entry. Please select 1, 2, 3, or 4.")
      }
    }
  }
  
  # Transform eFlow output to single row per sample/taxa
  # Identify the first column that matches a ps_sample_id value
  ps_sample_id_cols <- which(colnames(edna_data) %in% metadata_df$ps_sample_id)
  if (length(ps_sample_id_cols) > 0) {
    first_sample_col <- min(ps_sample_id_cols)
    id_vars <- colnames(edna_data)[1:(first_sample_col - 1)]
    
    # Reshape edna_data using melt from reshape2
    edna_data <- melt(edna_data, id.vars = id_vars, variable.name = "ps_sample_id", value.name = "sequence_count")
  } else {
    stop("No matching ps_sample_id columns found in edna_data")
  }
  
  
  
  # Merge metadata, labwork, and edna_data
  merged_data <- metadata_df %>% 
    left_join(labwork_df, by = "ps_sample_id") %>% 
    left_join(edna_data, by = "ps_sample_id")
  
  # Read FASTA file (DNA sequences) 
  # Revisit later to write out fasta files with identifier names Kate!
	if (!is.null(seq)) {
 	 dna_sequences <- readDNAStringSet(seq)
	}
  
  # Map columns to the occurrence data output format
  occurrence_df <- merged_data %>% 
    transmute(
      occurrenceID = ps_sample_id,
      source_mat_id = paste0(ps_sample_id, "_", OTU),
      modified = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
      basisOfRecord = "MaterialSample",
      recordedBy = team_lead,
      institutionID = "Pristine Seas",
      institutionCode = "",
      datasetName = expedition,
      parentEventID = expedition,
      eventID = paste(expedition, ps_station_id, depth_m, date, collection_time, sep = "_"),
      organismQuantityType = "DNA sequence reads",
      organismQuantity = "",
      occurrenceStatus = "present",
      associatedSequences = "",
      fieldNumber = ps_station_id,
      minimumDepthInMeters = depth_m,
      maximumDepthInMeters = depth_m,
      verbatimDepth = depth_m,
      eventDate = ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", date),date, format(as.Date(date, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%d-%b-%Y", "%b %d, %Y")),"%Y-%m-%d")),
      eventTime = collection_time,
      year = format(as.Date(date), "%Y"),
      month = format(as.Date(date), "%m"),
      day = format(as.Date(date), "%d"),
      verbatimEventDate = date,
      sampleSizeValue = water_liters,
      sampleSizeUnit = "liters",
      country = country,
      waterBody = waterBody,
      islandGroup = islandGroup,
      island = island,
      countryCode = countrycode::countrycode(country, "country.name", "iso2c"),
      locality = Atoll,
      decimalLatitude = lat,
      decimalLongitude = long,
      geodeticDatum = "WGS84",
      verbatimLatitude = lat,
      verbatimLongitude = long,
      preparations = paste(
        ifelse(!is.na(`Extraction Date`), paste("Extraction Date:", `Extraction Date`), NA),
        ifelse(!is.na(`Lysate (pre- extraction mL)`), paste("Lysate (pre- extraction mL):", `Lysate (pre- extraction mL)`), NA),
        ifelse(!is.na(`Extracted by`), paste("Extracted by:", `Extracted by`), NA),
        ifelse(!is.na(`Extraction Notes`), paste("Extraction Notes:", `Extraction Notes`), NA),
        ifelse(!is.na(`Extracted Qubit (ng/ul)`), paste("Extracted Qubit (ng/ul):", `Extracted Qubit (ng/ul)`), NA),
        ifelse(!is.na(`Extraction Amounts (uL)`), paste("Extraction Amounts (uL):", `Extraction Amounts (uL)`), NA),
        sep = " | "
      ),
      samplingProtocol = paste(
        ifelse(!is.na(preservative), paste("preservative:", preservative), NA),
        ifelse(!is.na(filter_type), paste("filter_type:", filter_type), NA),
        ifelse(!is.na(filter_time), paste("filter_time:", filter_time), NA),
        sep = " | "
      ),
      verbatimIdentification = "",
      scientificNameID = "",
      scientificName = ifelse(!is.na(genus) & !is.na(species), paste(genus, species),
                              ifelse(!is.na(species), species,
                                     coalesce(kingdom, phylum, class, order, family, genus, species))),
      higherClassification = domain,
      kingdom = kingdom,
      phylum = phylum,
      class = class,
      order = order,
      family = family,
      genus = genus,
      specificEpithet = species,
      taxonRank = names(which.max(!is.na(c(domain, kingdom, phylum, class, order, family, genus, specificEpithet)))),
      taxonID = taxid
    )
  
  # Map columns to the DNA-derived data output format
  dnaderiveddata_df <- merged_data %>% 
    transmute(
      source_mat_id = paste0(ps_sample_id, "_", OTU),
      project_name = expedition,
      env_broad_scale = env_broad_scale,
      env_medium = env_medium,
      lib_layout = lib_layout,
      target_gene = `Amplicon Target`,
      target_subfragment = "",
      seq_meth = seq_meth,
      pcr_cond = paste(
        ifelse(!is.na(`PCR Date`), paste("PCR Date:", `PCR Date`), NA),
        ifelse(!is.na(`DNA template Quantity (ul)`), paste("DNA template Quantity (ul):", `DNA template Quantity (ul)`), NA),
        ifelse(!is.na(`Primer Quantity (ul)`), paste("Primer Quantity (ul):", `Primer Quantity (ul)`), NA),
        ifelse(!is.na(`PCR Notes`), paste("PCR Notes:", `PCR Notes`), NA),
        ifelse(!is.na(`Post Bead Clean Qubit (ng/ul)`), paste("Post Bead Clean Qubit (ng/ul):", `Post Bead Clean Qubit (ng/ul)`), NA),
        ifelse(!is.na(`Adaptor Tag`), paste("Adaptor Tag:", `Adaptor Tag`), NA),
        ifelse(!is.na(`Post Library Qubit (ng/ul)`), paste("Post Library Qubit (ng/ul):", `Post Library Qubit (ng/ul)`), NA),
        ifelse(!is.na(`BioAnalzer / qPCR`), paste("BioAnalzer / qPCR:", `BioAnalzer / qPCR`), NA),
        sep = " | "
      ),
      pcr_primers = paste0("FWD:",`F-Primer`,";REV:",`R-Primer`),
      DNA_sequence = sequence
    )
  
  
  # Create occurrence.txt and dnaderiveddata.txt files
  write_tsv(occurrence_df, "occurrence.txt", na="")
  write_tsv(dnaderiveddata_df, "dnaderiveddata.txt", na="")
  
  # Generate eml.xml with purpose, intellectual rights, and associated party
  eml_doc <- newXMLNode("eml:eml", namespaceDefinitions = c(
    "eml" = "eml://ecoinformatics.org/eml-2.1.1"
  )
  )
  dataset <- newXMLNode("dataset", parent = eml_doc)
  newXMLNode("purpose", "These data were made accessible through Pristine Seas.", parent = dataset)
  intellectual_rights <- newXMLNode("intellectualRights", parent = dataset)
  newXMLNode("para", "To the extent possible under law, the publisher has waived all rights to these data and has dedicated them to the", parent = intellectual_rights)
  newXMLNode("ulink", attrs = c(url = "http://creativecommons.org/publicdomain/zero/1.0/legalcode"), "Public Domain (CC0 1.0)", parent = intellectual_rights)
  
  # Add associated party
  associated_party <- newXMLNode("associatedParty", parent = dataset)
  newXMLNode("individualName", metadata_df$team_lead[1], parent = associated_party)
  newXMLNode("organizationName", "Pristine Seas", parent = associated_party)
  
  saveXML(eml_doc, file = "eml.xml")
  
  # Generate meta.xml with dynamic field indexing using XML package
  meta_doc <- newXMLNode("archive", 
                         namespaceDefinitions = c(
                           "dwc" = "http://rs.tdwg.org/dwc/text/"
                         ),
                         attrs = c("metadata" = "eml.xml")
  )
  
  # Core section (Occurrence data)
  core <- newXMLNode("core", 
                     parent = meta_doc,
                     attrs = c(
                       "encoding" = "UTF-8", 
                       "fieldsTerminatedBy" = "\t", 
                       "linesTerminatedBy" = "\n", 
                       "fieldsEnclosedBy" = "", 
                       "ignoreHeaderLines" = "1", 
                       "rowType" = "http://rs.tdwg.org/dwc/terms/Occurrence"
                     )
  )
  newXMLNode("files", newXMLNode("location", "occurrence.txt"), parent = core)
  newXMLNode("id", attrs = c("index" = "0"), parent = core)
  
  for (i in seq_along(colnames(occurrence_df))) {
    newXMLNode("field", 
               attrs = c("index" = as.character(i), "term" = paste0("http://rs.tdwg.org/dwc/terms/", colnames(occurrence_df)[i])),
               parent = core
    )
  }
  
  # Extension section (DNA-derived data)
  extension <- newXMLNode("extension", 
                          parent = meta_doc,
                          attrs = c(
                            "encoding" = "UTF-8", 
                            "fieldsTerminatedBy" = "\t", 
                            "linesTerminatedBy" = "\n", 
                            "fieldsEnclosedBy" = "", 
                            "ignoreHeaderLines" = "1", 
                            "rowType" = "http://rs.gbif.org/terms/1.0/DNADerivedData"
                          )
  )
  newXMLNode("files", newXMLNode("location", "dnaderiveddata.txt"), parent = extension)
  newXMLNode("coreid", attrs = c("index" = "0"), parent = extension)
  
  mixs_terms <- list(
    "source_mat_id" = "0000026",
    "project_name" = "0000092",
    "env_broad_scale" = "0000012",
    "env_medium" = "0000014",
    "lib_layout" = "0000041",
    "target_gene" = "0000044",
    "target_subfragment" = "0000045",
    "seq_meth" = "0000050",
    "pcr_cond" = "0000049",
    "pcr_primers" = "0000046"
  )
  
  for (i in seq_along(colnames(dnaderiveddata_df))) {
    col_name <- colnames(dnaderiveddata_df)[i]
    if (col_name != "DNA_sequence") {
      term_url <- ifelse(col_name %in% names(mixs_terms), 
                         paste0("https://genomicsstandardsconsortium.github.io/mixs/", mixs_terms[[col_name]]), 
                         paste0("http://rs.gbif.org/terms/", col_name)
      )
      newXMLNode("field", 
                 attrs = c("index" = as.character(i), "term" = term_url),
                 parent = extension
      )
    }
  }
  
  # Save meta.xml file
  saveXML(meta_doc, file = "meta.xml")
  
  # Define the output file name in the working directory
  output_zip <- file.path(getwd(), "dwca_output.zip")
  
  # Ensure all necessary files exist before zipping
  files_to_zip <- c("occurrence.txt", "dnaderiveddata.txt", "meta.xml", "eml.xml")
  files_to_zip <- files_to_zip[file.exists(files_to_zip)]
  
  # Zip the files using zipr (ensures correct writing to the working directory)
  zip::zipr(zipfile = output_zip, files = files_to_zip)
  
  # Confirm the zip file is created
  if (file.exists(output_zip)) {
    message("DwC-a file successfully created at: ", output_zip)
  } else {
    stop("Error: Failed to create DwC-a file")
  }
}
