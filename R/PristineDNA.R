#Check and install required packages
required_packages <- c("magrittr","countrycode", "measurements", "readxl", "readr", "dplyr", 
                       "stringr", "Biostrings", "reshape2", "purrr", "XML", "utils", "zip",
                       "lubridate")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages, dependencies = TRUE)
}

# Load  libraries
library(magrittr) #fixing a dplyr issue
library(countrycode)
library(measurements)
library(readxl)
library(readr) #may not be necessary anymore - KS check
library(dplyr)
library(stringr)
library(Biostrings)
library(reshape2)
library(purrr)
library(XML)
library(utils)
library(zip)
library(lubridate)

#############################

#function to convert latitude and longitude to decimal degrees
convert_to_decimal <- function(coord) {
  if (is.na(coord) | coord == "") return(NA)
  
  coord <- gsub("[°]", " ", coord)
  coord <- gsub("[’']", " ", coord)
  coord <- gsub('["]', " ", coord)
  
  parts <- unlist(strsplit(trimws(coord), "\\s+"))
  direction <- parts[length(parts)]
  numeric_parts <- as.numeric(parts[1:(length(parts) - 1)])
  
  decimal <- numeric_parts[1] + 
    (ifelse(length(numeric_parts) > 1, numeric_parts[2] / 60, 0)) +
    (ifelse(length(numeric_parts) > 2, numeric_parts[3] / 3600, 0))
  
  if (direction %in% c("S", "W")) {
    decimal <- -decimal
  }
  
  return(decimal)
}

###############################

#function to parse and format any date input into ISO8601
format_to_iso8601 <- function(date_input) {
  parsed_date <- as.Date(lubridate::parse_date_time(date_input, orders = c("ymd", "mdy", "dmy", "ydm")))
  return(format(parsed_date, "%Y-%m-%d"))
}

################################

######placeholder to add 24h clock time conversion

#########################################################################################################

#Major function to load eDNA data and map to occurrence and DNA-derived data output format
#Basic steps are reading in standard lab files
#performing initial QC filtering based on OBIS/GBIF minimal data requirements

pristinedna <- function(metadata, eflow, seq = NULL) {
  
  #read in Excel file (metadata and labwork sheets)
  metadata_df <- read_excel(metadata, sheet = "metadata")
  labwork_df <- read_excel(metadata, sheet = "labwork")
  
  #read in TSV file (eFlow output) -- KS CONFIRM WITH LAB IF THIS IS AUTOMATED
  edna_data <- read.delim(eflow, na = c("NA", "dropped"),header=TRUE)
  
  ######Start pre-processing QC and conversion steps
  
  #QC check for missing ps_sample_id values as this is a fatal OBIS/GBIF error later and confuses data melt/merging
  missing_metadata <- metadata_df %>% filter(is.na(ps_sample_id) | ps_sample_id == "" | ps_sample_id == 0)
  missing_labwork <- labwork_df %>% filter(is.na(ps_sample_id) | ps_sample_id == "" | ps_sample_id == 0)
  
  if (nrow(missing_metadata) > 0) {
    message(paste(nrow(missing_metadata), "rows have missing ps_sample_id values in the 'metadata' sheet."))
  }
  if (nrow(missing_labwork) > 0) {
    message(paste(nrow(missing_labwork), "rows have missing ps_sample_id values in the 'labwork' sheet."))
  }
  
  missing_count <- nrow(missing_metadata) + nrow(missing_labwork)
  if (missing_count > 0) {
    repeat {
      user_choice <- readline(paste("Options:\n",
                                    "1 - Review rows now and decide what to do\n",
                                    "2 - Remove rows without reviewing them\n",
                                    "3 - Quit program and manually remove these\n",
                                    "Enter option: "))
      if (user_choice %in% c("1", "2", "3")) break
      message("Invalid option. Please enter 1, 2, or 3.")
    }
    
    if (user_choice == "3") {
      stop("Good luck!")
    } else if (user_choice == "2") {
      metadata_df <- metadata_df %>% filter(!is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
      labwork_df <- labwork_df %>% filter(!is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
    } else if (user_choice == "1") {
      if (nrow(missing_metadata) > 0) {
        print("Reviewing missing ps_sample_id rows from 'metadata', 10 at a time:")
        for (i in seq(1, nrow(missing_metadata), by = 10)) {
          print(missing_metadata[i:min(i+9, nrow(missing_metadata)), ])
          if (i + 9 < nrow(missing_metadata)) {
            cont <- readline("Press Enter to see more, or type 'exit' to stop reviewing: ")
            if (tolower(cont) == "exit") break
          }
        }
      }
      if (nrow(missing_labwork) > 0) {
        print("Reviewing missing ps_sample_id rows from 'labwork', 10 at a time:")
        for (i in seq(1, nrow(missing_labwork), by = 10)) {
          print(missing_labwork[i:min(i+9, nrow(missing_labwork)), ])
          if (i + 9 < nrow(missing_labwork)) {
            cont <- readline("Press Enter to see more, or type 'exit' to stop reviewing: ")
            if (tolower(cont) == "exit") break
          }
        }
      }
      repeat {
        final_choice <- readline("Options:\n1 - OK to remove rows from ongoing formatting\n2 - Quit program and manually revise\nEnter option: ")
        if (final_choice %in% c("1", "2")) break
        message("Invalid option. Please enter 1 or 2.")
      }
      if (final_choice == "2") stop("Good luck!")
      metadata_df <- metadata_df %>% filter(!is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
      labwork_df <- labwork_df %>% filter(!is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
    }
  }
  
  #convert lat/long to decimal degrees
  metadata_df$lat <- sapply(metadata_df$lat, convert_to_decimal)
  metadata_df$long <- sapply(metadata_df$long, convert_to_decimal)
  
  #convert dates to ISO 8601
  metadata_df$date <- sapply(metadata_df$date , format_to_iso8601)
  labwork_df$`Extraction Date` <- sapply(labwork_df$`Extraction Date`, format_to_iso8601)
  
  #convert - hold for entry to convert time to 24h clock
  #Needs to happen for metadata_df$collection_time and labwork_df$filter_time
  
  #QC Check for pid values <90 and give user choice to drop
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
  
  #QC Check for complete missing taxonomic information to drop
  #empty/dropped taxonomy will result in empty scientificName which is a fatal OBIS/GBIF error
  
  taxonomic_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  #find, store and remove rows where all taxonomic columns are missing or contain unwanted values
  missing_taxa_rows <- edna_data %>% 
    filter(dplyr::if_all(all_of(taxonomic_cols), ~ is.na(.) | . == "" | . == "0" | . == "dropped"))
  
  if (nrow(missing_taxa_rows) > 0) {
    removed_rows <- missing_taxa_rows  # Save before filtering
  }
  edna_data <- edna_data %>% 
    filter(!(dplyr::if_all(all_of(taxonomic_cols), ~ is.na(.) | . == "" | . == "0" | . == "dropped")))
  
  #print removed rows after filtering to screen and save as CSV file 
  if (exists("removed_rows") && nrow(removed_rows) > 0) {
    message(paste(nrow(removed_rows), "rows have been removed due to missing taxonomic information. This has been written to removed_eFlow_rows.csv"))
    write.csv(removed_rows, file = "removed_eFlow_rows.csv", row.names = FALSE)
  }
  
  #reshape edna_data to get it ready for melting of only the columns we want
  #identify the first occurrence of a ps_sample_id column, there may be missing or extra samples
  ps_sample_id_cols <- which(colnames(edna_data) %in% metadata_df$ps_sample_id)
  
  #preserve all columns in the eFlow/eDNA file that are to the left of the first sample count instance.
  #example files include PCR values to the right which we aren't using here.
  #melt data for all metadata columns to the left of the count data, and only count data that has shared ps_sample_id
  
  if (length(ps_sample_id_cols) > 0) {
    first_sample_col <- min(ps_sample_id_cols)
    id_vars <- colnames(edna_data)[1:(first_sample_col - 1)]
    
    edna_data <- melt(edna_data, id.vars = id_vars, measure.vars = colnames(edna_data)[ps_sample_id_cols], 
                      variable.name = "ps_sample_id", value.name = "sequence_count")
  } else {
    stop("No matching ps_sample_id columns found in edna_data")
  }
  
  #confirm ps_sample_id is in output as this is a fatal error
  if (!"ps_sample_id" %in% colnames(edna_data)) {
    stop("Error: ps_sample_id column is missing after melting.")
  }
  
  #check that all 3 datasets have the same ps_sample_ids before merge. If not, user decision.
  common_ps_sample_ids <- Reduce(intersect, list(metadata_df$ps_sample_id, labwork_df$ps_sample_id, edna_data$ps_sample_id))
  
  #identify and print to user
  missing_metadata <- setdiff(metadata_df$ps_sample_id, common_ps_sample_ids)
  missing_labwork <- setdiff(labwork_df$ps_sample_id, common_ps_sample_ids)
  missing_edna <- setdiff(edna_data$ps_sample_id, common_ps_sample_ids)
  
  missing_count_metadata <- length(missing_metadata)
  missing_count_labwork <- length(missing_labwork)
  missing_count_edna <- length(missing_edna)
  
  total_missing <- missing_count_metadata + missing_count_labwork + missing_count_edna
  
  if (total_missing > 0) {
    message(paste(missing_count_metadata, "missing samples present in metadata."))
    message(paste(missing_count_labwork, "missing samples present in labwork."))
    message(paste(missing_count_edna, "missing samples present in eDNA data."))
    
    repeat {
      user_choice <- readline(paste("Options:\n",
                                    "1 - Remove samples missing from any dataset and continue\n",
                                    "2 - Review missing samples\n",
                                    "3 - Quit program and manually fix\n",
                                    "Enter option: "))
      if (user_choice %in% c("1", "2", "3")) break
      message("Invalid option. Please enter 1, 2, or 3.")
    }
    
    if (user_choice == "3") stop("Good luck!")
    
    if (user_choice == "2") {
      message("Reviewing missing samples by dataset:")
      if (missing_count_metadata > 0) {
        print("Missing samples present in metadata:")
        print(metadata_df %>% filter(ps_sample_id %in% missing_metadata))
      } else if (missing_count_edna > 0) {
        print("Missing samples present in eDNA data:")
        print(edna_data %>% filter(ps_sample_id %in% missing_edna))
      } else if (missing_count_labwork > 0) {
        print("Missing samples present in labwork:")
        print(labwork_df %>% filter(ps_sample_id %in% missing_labwork))
      }
      repeat {
        final_choice <- readline("Options:\n1 - OK to remove missing samples and continue\n2 - Quit program and manually revise\nEnter option: ")
        if (final_choice %in% c("1", "2")) break
        message("Invalid option. Please enter 1 or 2.")
      }
      if (final_choice == "2") stop("Good luck!")
    }
    
    #remove samples not present in all three datasets
    metadata_df <- metadata_df %>% filter(ps_sample_id %in% common_ps_sample_ids)
    labwork_df <- labwork_df %>% filter(ps_sample_id %in% common_ps_sample_ids)
    edna_data <- edna_data %>% filter(ps_sample_id %in% common_ps_sample_ids)
  }
  
  #perform inner join to ensure only common ps_sample_id values are merged
  merged_data <- metadata_df %>% 
    dplyr::inner_join(labwork_df, by = "ps_sample_id") %>% 
    dplyr::inner_join(edna_data, by = "ps_sample_id")
  
  
  
  # Read FASTA file (DNA sequences) 
  # Revisit later to write out fasta files with identifier names Kate!
  if (!is.null(seq)) {
    dna_sequences <- readDNAStringSet(seq)
  }
  
  ######Convert and map columns to OBIS/GBIF-approved values
  
  #occurrence file output
  occurrence_df <- merged_data %>% 
    transmute(
      occurrenceID = paste(ps_sample_id, OTU, sep = "_"),
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
      organismQuantity = sequence_count,
      occurrenceStatus = ifelse(sequence_count > 0, "present", "absent"),
      associatedSequences = "",
      fieldNumber = ps_station_id,
      minimumDepthInMeters = depth_m,
      maximumDepthInMeters = depth_m,
      verbatimDepth = depth_m,
      eventDate = format(as.Date(date), "%Y-%m-%d"),
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
                                     coalesce(domain, kingdom, phylum, class, order, family, genus, species))),
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
  
  #dnaderived data file output
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
  
  
  #save occurrence.txt and dnaderiveddata.txt to local directory
  write_tsv(occurrence_df, "occurrence.txt", na="")
  write_tsv(dnaderiveddata_df, "dnaderiveddata.txt", na="")
  
  #generate eml.xml with purpose, intellectual rights, and associated party
  eml_doc <- newXMLNode("eml:eml", namespaceDefinitions = c(
    "eml" = "eml://ecoinformatics.org/eml-2.1.1"
  )
  )
  dataset <- newXMLNode("dataset", parent = eml_doc)
  newXMLNode("purpose", "These data were made accessible through Pristine Seas.", parent = dataset)
  intellectual_rights <- newXMLNode("intellectualRights", parent = dataset)
  newXMLNode("para", "To the extent possible under law, the publisher has waived all rights to these data and has dedicated them to the", parent = intellectual_rights)
  newXMLNode("ulink", attrs = c(url = "http://creativecommons.org/publicdomain/zero/1.0/legalcode"), "Public Domain (CC0 1.0)", parent = intellectual_rights)
  
  #associated party
  associated_party <- newXMLNode("associatedParty", parent = dataset)
  newXMLNode("individualName", metadata_df$team_lead[1], parent = associated_party)
  newXMLNode("organizationName", "Pristine Seas", parent = associated_party)
  
  saveXML(eml_doc, file = "eml.xml")
  
  #generate meta.xml with dynamic field indexing using XML package
  meta_doc <- newXMLNode("archive", 
                         namespaceDefinitions = c(
                           "dwc" = "http://rs.tdwg.org/dwc/text/"
                         ),
                         attrs = c("metadata" = "eml.xml")
  )
  
  #core section for occurrence data
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
  
  #dnaderived data extension section
  #note this is a relatively new extension and terms are not well used or defined. References currently have to be hardcoded!
  #revisit these terms in future to see if code can be generalized
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
  
  #save meta.xml file
  saveXML(meta_doc, file = "meta.xml")
  
  #define the output file name in the working directory
  output_zip <- file.path(getwd(), "dwca_output.zip")
  
  #ensure all necessary files exist then zip
  files_to_zip <- c("occurrence.txt", "dnaderiveddata.txt", "meta.xml", "eml.xml")
  files_to_zip <- files_to_zip[file.exists(files_to_zip)]
  zip::zipr(zipfile = output_zip, files = files_to_zip)
  
  #confirm zip file is created and save to local directory
  if (file.exists(output_zip)) {
    message(paste("DwC-A file successfully created at", output_zip, 
                  "- Please review your files to ensure they look correct, and run them through pristineQC before submitting to OBIS/GBIF."))
  } else {
    stop("Error: Failed to create DwC-A file")
  }
}
