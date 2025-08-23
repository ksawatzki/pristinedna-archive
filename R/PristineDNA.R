# Set up
# function to convert latitude and longitude to decimal degrees flexibly
convert_to_decimal <- function(x) {
  vapply(x, FUN.VALUE = numeric(1), function(coord) {
    if (is.null(coord) || is.na(coord) || coord == "") return(NA_real_)
    coord <- trimws(as.character(coord))
    dec_try <- suppressWarnings(as.numeric(coord))
    if (!is.na(dec_try)) return(dec_try)
    up <- toupper(coord)
    dir_match <- regmatches(up, regexpr("[NSEW](?!.*[NSEW])", up, perl = TRUE))
    has_dir <- length(dir_match) == 1
    dir <- if (has_dir) dir_match else NA_character_
    norm <- coord
    norm <- gsub("[°º∘]", " ", norm)     # degrees
    norm <- gsub("[′’ʹ']", " ", norm)    # minutes
    norm <- gsub("[″\"ʺ]", " ", norm)    # seconds
    norm <- gsub(",", " ", norm)
    norm <- gsub("[A-Za-z]", " ", norm)
    norm <- gsub("\\s+", " ", norm)
    parts <- strsplit(trimws(norm), " ", fixed = FALSE)[[1]]
    nums  <- suppressWarnings(as.numeric(parts))
    nums  <- nums[!is.na(nums)]
    if (length(nums) == 0) return(NA_real_)
    deg <- nums[1]
    min <- if (length(nums) >= 2) nums[2] else 0
    sec <- if (length(nums) >= 3) nums[3] else 0
    dec <- abs(deg) + min/60 + sec/3600
    sgn <- if (deg < 0) -1 else 1
    if (has_dir) {
      if (dir %in% c("S", "W")) sgn <- -1
      if (dir %in% c("N", "E")) sgn <-  1
    }
    sgn * dec
  })
}
# function to parse and format any date input into ISO8601
format_to_iso8601 <- function(date_input) {
  parsed_date <- as.Date(lubridate::parse_date_time(date_input, orders = c("ymd", "mdy", "dmy", "ydm")))
  return(format(parsed_date, "%Y-%m-%d"))
}

#' Process eDNA data for DwC-A submission (occurrence + DNA-derived data)
#'
#' Reads your Pristine Seas metadata workbook (`metadata`), eFlow TSV (`eflow`),
#' and (optionally) a FASTA file of sequences, applies QC checks, and writes a
#' Darwin Core Archive (occurrence.txt, dnaderiveddata.txt, meta.xml, eml.xml)
#' plus a `dwca_output.zip` in the working directory.
#'
#' Interactive prompts help resolve missing IDs, low % identity rows, and sample
#' mismatches across inputs when `qc=TRUE`.
#'
#' @param metadata Path to the Excel workbook containing sheets `metadata` and `labwork`.
#' @param eflow Path to eFlow TSV (wide counts table + taxonomy columns).
#' @param seq Optional path to FASTA file of DNA sequences.
#' @param qc Logical; if `TRUE`, run OBIS + Hmisc QC and write summaries.
#'
#' @return Invisibly returns a list with `occurrence_df` and `dnaderiveddata_df`.
#'
#' @examples
#' \dontrun{
#' pristinedna("my_metadata.xlsx", "eflow_output.tsv", qc = TRUE)
#' }
#'
#' @importFrom readxl read_excel
#' @importFrom readr read_tsv write_tsv
#' @importFrom dplyr filter mutate left_join inner_join transmute
#' @importFrom stringr str_detect str_replace str_trim
#' @importFrom reshape2 melt
#' @importFrom Biostrings readDNAStringSet
#' @importFrom XML newXMLNode saveXML
#' @importFrom zip zipr
#' @importFrom countrycode countrycode
#' @importFrom utils read.delim
#' @importFrom magrittr %>%
#' @export
pristinedna <- function(metadata, eflow, seq = NULL, qc = TRUE) {
  # read in Excel file (metadata and labwork sheets)
  metadata_df <- readxl::read_excel(metadata, sheet = "metadata")
  labwork_df <- readxl::read_excel(metadata, sheet = "labwork")
  
  # read in TSV file (eFlow output)
  edna_data <- utils::read.delim(eflow, na = c("NA", "dropped"), header = TRUE)
  
  ###### Start pre-processing QC and conversion steps
  
  # QC check for missing ps_sample_id values
  missing_metadata <- dplyr::filter(metadata_df, is.na(ps_sample_id) | ps_sample_id == "" | ps_sample_id == 0)
  missing_labwork <- dplyr::filter(labwork_df, is.na(ps_sample_id) | ps_sample_id == "" | ps_sample_id == 0)
  
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
      metadata_df <- dplyr::filter(metadata_df, !is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
      labwork_df <- dplyr::filter(labwork_df, !is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
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
      metadata_df <- dplyr::filter(metadata_df, !is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
      labwork_df <- dplyr::filter(labwork_df, !is.na(ps_sample_id) & ps_sample_id != "" & ps_sample_id != 0)
    }
  }
  
  # convert lat/long to decimal degrees
  metadata_df$lat  <- sapply(metadata_df$lat,  convert_to_decimal)
  metadata_df$long <- sapply(metadata_df$long, convert_to_decimal)
  
  # convert dates to ISO 8601
  metadata_df$date               <- sapply(metadata_df$date, format_to_iso8601)
  labwork_df$`Extraction Date`   <- sapply(labwork_df$`Extraction Date`, format_to_iso8601)
  
  # QC Check for pid values <90 and user choice to drop
  low_pid_rows <- dplyr::filter(edna_data, pid < 90)
  
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
        edna_data <- dplyr::filter(edna_data, pid >= 90)
        break
      } else if (user_choice == "4") {
        stop("See you later!")
      } else if (user_choice == "3") {
        for (i in 1:nrow(low_pid_rows)) {
          print(low_pid_rows[i, c("genus", "species", "sequence", "seq_length", "OTU", "unique_hits", "pid")])
          repeat {
            decision <- readline("1 - Keep\n2 - Drop\n3 - Quit program\nEnter choice: ")
            if (decision == "1") {
              break
            } else if (decision == "3") {
              stop("See you later!")
            } else if (decision == "2") {
              edna_data <- dplyr::filter(edna_data, !(pid == low_pid_rows$pid[i] & genus == low_pid_rows$genus[i] & species == low_pid_rows$species[i]))
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
  
  # QC for missing taxonomy (fatal for OBIS/GBIF)
  taxonomic_cols <- c("domain","kingdom","phylum","class","order","family","genus","species")
  missing_taxa_rows <- dplyr::filter(edna_data, dplyr::if_all(dplyr::all_of(taxonomic_cols), ~ is.na(.) | . == "" | . == "0" | . == "dropped"))
  if (nrow(missing_taxa_rows) > 0) {
    removed_rows <- missing_taxa_rows
  }
  edna_data <- dplyr::filter(edna_data, !(dplyr::if_all(dplyr::all_of(taxonomic_cols), ~ is.na(.) | . == "" | . == "0" | . == "dropped")))
  if (exists("removed_rows") && nrow(removed_rows) > 0) {
    message(paste(nrow(removed_rows), "rows have been removed due to missing taxonomic information. This has been written to removed_eFlow_rows.csv"))
    utils::write.csv(removed_rows, file = "removed_eFlow_rows.csv", row.names = FALSE)
  }
  
  # Reshape eFlow: melt counts for ps_sample_id columns present in metadata
  ps_sample_id_cols <- which(colnames(edna_data) %in% metadata_df$ps_sample_id)
  if (length(ps_sample_id_cols) > 0) {
    first_sample_col <- min(ps_sample_id_cols)
    id_vars <- colnames(edna_data)[1:(first_sample_col - 1)]
    edna_data <- reshape2::melt(edna_data, id.vars = id_vars, measure.vars = colnames(edna_data)[ps_sample_id_cols],
                                variable.name = "ps_sample_id", value.name = "sequence_count")
  } else {
    stop("No matching ps_sample_id columns found in edna_data")
  }
  
  if (!"ps_sample_id" %in% colnames(edna_data)) {
    stop("Error: ps_sample_id column is missing after melting.")
  }
  
  # Enforce intersection across three inputs
  common_ps_sample_ids <- Reduce(intersect, list(metadata_df$ps_sample_id, labwork_df$ps_sample_id, edna_data$ps_sample_id))
  missing_metadata <- setdiff(metadata_df$ps_sample_id, common_ps_sample_ids)
  missing_labwork  <- setdiff(labwork_df$ps_sample_id,  common_ps_sample_ids)
  missing_edna     <- setdiff(edna_data$ps_sample_id,   common_ps_sample_ids)
  total_missing    <- length(missing_metadata) + length(missing_labwork) + length(missing_edna)
  
  if (total_missing > 0) {
    message(paste(length(missing_metadata), "missing samples present in metadata."))
    message(paste(length(missing_labwork),  "missing samples present in labwork."))
    message(paste(length(missing_edna),     "missing samples present in eDNA data."))
    
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
      if (length(missing_metadata) > 0) {
        print("Missing samples present in metadata:"); print(dplyr::filter(metadata_df, ps_sample_id %in% missing_metadata))
      } else if (length(missing_edna) > 0) {
        print("Missing samples present in eDNA data:"); print(dplyr::filter(edna_data, ps_sample_id %in% missing_edna))
      } else if (length(missing_labwork) > 0) {
        print("Missing samples present in labwork:"); print(dplyr::filter(labwork_df, ps_sample_id %in% missing_labwork))
      }
      repeat {
        final_choice <- readline("Options:\n1 - OK to remove missing samples and continue\n2 - Quit program and manually revise\nEnter option: ")
        if (final_choice %in% c("1", "2")) break
        message("Invalid option. Please enter 1 or 2.")
      }
      if (final_choice == "2") stop("Good luck!")
    }
    
    metadata_df <- dplyr::filter(metadata_df, ps_sample_id %in% common_ps_sample_ids)
    labwork_df  <- dplyr::filter(labwork_df,  ps_sample_id %in% common_ps_sample_ids)
    edna_data   <- dplyr::filter(edna_data,   ps_sample_id %in% common_ps_sample_ids)
  }
  
  # Merge
  merged_data <- dplyr::inner_join(metadata_df, labwork_df, by = "ps_sample_id") %>%
    dplyr::inner_join(edna_data,   by = "ps_sample_id")
  
  # Optional FASTA
  if (!is.null(seq)) {
    dna_sequences <- Biostrings::readDNAStringSet(seq)
  }
  
  ###### Map to DwC terms
  
  occurrence_df <- dplyr::transmute(
    merged_data,
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
                                   dplyr::coalesce(domain, kingdom, phylum, class, order, family, genus, species))),
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
  
  dnaderiveddata_df <- dplyr::transmute(
    merged_data,
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
  
  # QC via obistools/Hmisc
  if (qc) {
    message("\n=== Running QC checks with obistools and Hmisc ===")
    if (!requireNamespace("obistools", quietly = TRUE)) stop("Package 'obistools' is required for QC but not installed.")
    if (!requireNamespace("Hmisc", quietly = TRUE))     stop("Package 'Hmisc' is required for QC but not installed.")
    suppressWarnings(library(obistools))
    suppressWarnings(library(Hmisc))
    
    message("\n-- OBIS Field Checks --")
    fields_check <- obistools::check_fields(occurrence_df)
    event_check  <- tryCatch(obistools:::check_eventdate(occurrence_df), error = function(e) data.frame())
    coord_check  <- tryCatch(obistools:::check_coordinates(occurrence_df), error = function(e) data.frame())
    
    message(paste("-", nrow(fields_check), "rows with invalid DwC fields"))
    message(paste("-", nrow(event_check),  "rows with problematic eventDate"))
    message(paste("-", nrow(coord_check),  "rows with invalid coordinates"))
    
    if (nrow(coord_check) > 0) {
      utils::write.csv(coord_check, "bad_coordinates.csv", row.names = FALSE)
      message("→ Saved bad coordinates to 'bad_coordinates.csv'")
    }
    if (nrow(event_check) > 0) {
      utils::write.csv(event_check, "bad_event_dates.csv", row.names = FALSE)
      message("→ Saved problematic eventDate entries to 'bad_event_dates.csv'")
    }
    if (nrow(fields_check) > 0) {
      utils::write.csv(fields_check, "bad_fields.csv", row.names = FALSE)
      message("→ Saved invalid DwC field rows to 'bad_fields.csv'")
    }
    
    message("\n-- Hmisc Variable Summaries --")
    sink("qc_summary_occurrence.txt"); print(Hmisc::describe(occurrence_df)); sink()
    message("→ Variable summary for occurrence_df written to 'qc_summary_occurrence.txt'")
    sink("qc_summary_dnaderiveddata.txt"); print(Hmisc::describe(dnaderiveddata_df)); sink()
    message("→ Variable summary for dnaderiveddata_df written to 'qc_summary_dnaderiveddata.txt'")
  }
  
  # Write outputs
  readr::write_tsv(occurrence_df, "occurrence.txt", na = "")
  readr::write_tsv(dnaderiveddata_df, "dnaderiveddata.txt", na = "")
  
  # Build EML
  eml_doc <- XML::newXMLNode("eml:eml", namespaceDefinitions = c("eml" = "eml://ecoinformatics.org/eml-2.1.1"))
  dataset <- XML::newXMLNode("dataset", parent = eml_doc)
  XML::newXMLNode("purpose", "These data were made accessible through Pristine Seas.", parent = dataset)
  intellectual_rights <- XML::newXMLNode("intellectualRights", parent = dataset)
  XML::newXMLNode("para", "To the extent possible under law, the publisher has waived all rights to these data and has dedicated them to the", parent = intellectual_rights)
  XML::newXMLNode("ulink", attrs = c(url = "http://creativecommons.org/publicdomain/zero/1.0/legalcode"), "Public Domain (CC0 1.0)", parent = intellectual_rights)
  associated_party <- XML::newXMLNode("associatedParty", parent = dataset)
  XML::newXMLNode("individualName", metadata_df$team_lead[1], parent = associated_party)
  XML::newXMLNode("organizationName", "Pristine Seas", parent = associated_party)
  XML::saveXML(eml_doc, file = "eml.xml")
  
  # Build meta.xml
  meta_doc <- XML::newXMLNode("archive",
                              namespaceDefinitions = c("dwc" = "http://rs.tdwg.org/dwc/text/"),
                              attrs = c("metadata" = "eml.xml"))
  core <- XML::newXMLNode("core", parent = meta_doc,
                          attrs = c("encoding"="UTF-8","fieldsTerminatedBy"="\t","linesTerminatedBy"="\n",
                                    "fieldsEnclosedBy"="","ignoreHeaderLines"="1",
                                    "rowType"="http://rs.tdwg.org/dwc/terms/Occurrence"))
  XML::newXMLNode("files", XML::newXMLNode("location", "occurrence.txt"), parent = core)
  XML::newXMLNode("id", attrs = c("index" = "0"), parent = core)
  for (i in seq_along(colnames(occurrence_df))) {
    XML::newXMLNode("field", attrs = c("index" = as.character(i),
                                       "term" = paste0("http://rs.tdwg.org/dwc/terms/", colnames(occurrence_df)[i])), parent = core)
  }
  
  extension <- XML::newXMLNode("extension", parent = meta_doc,
                               attrs = c("encoding"="UTF-8","fieldsTerminatedBy"="\t","linesTerminatedBy"="\n",
                                         "fieldsEnclosedBy"="","ignoreHeaderLines"="1",
                                         "rowType"="http://rs.gbif.org/terms/1.0/DNADerivedData"))
  XML::newXMLNode("files", XML::newXMLNode("location", "dnaderiveddata.txt"), parent = extension)
  XML::newXMLNode("coreid", attrs = c("index" = "0"), parent = extension)
  
  mixs_terms <- list("source_mat_id"="0000026","project_name"="0000092","env_broad_scale"="0000012",
                     "env_medium"="0000014","lib_layout"="0000041","target_gene"="0000044","target_subfragment"="0000045",
                     "seq_meth"="0000050","pcr_cond"="0000049","pcr_primers"="0000046")
  for (i in seq_along(colnames(dnaderiveddata_df))) {
    col_name <- colnames(dnaderiveddata_df)[i]
    if (col_name != "DNA_sequence") {
      term_url <- ifelse(col_name %in% names(mixs_terms),
                         paste0("https://genomicsstandardsconsortium.github.io/mixs/", mixs_terms[[col_name]]),
                         paste0("http://rs.gbif.org/terms/", col_name))
      XML::newXMLNode("field", attrs = c("index" = as.character(i), "term" = term_url), parent = extension)
    }
  }
  
  XML::saveXML(meta_doc, file = "meta.xml")
  
  # Zip
  output_zip <- file.path(getwd(), "dwca_output.zip")
  files_to_zip <- c("occurrence.txt","dnaderiveddata.txt","meta.xml","eml.xml")
  files_to_zip <- files_to_zip[file.exists(files_to_zip)]
  zip::zipr(zipfile = output_zip, files = files_to_zip)
  
  if (file.exists(output_zip)) {
    message(paste("DwC-A file successfully created at", output_zip))
  } else {
    stop("Error: Failed to create DwC-A file")
  }
  
  invisible(list(occurrence_df = occurrence_df, dnaderiveddata_df = dnaderiveddata_df))
}
