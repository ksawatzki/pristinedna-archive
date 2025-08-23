#' Generate SRA submission files from Pristine Seas metadata
#'
#' Interactively prepares two tab-delimited files for SRA submission:
#' a BioSample table tailored to MIMS/MIMAG/MIMARKS and an SRAmetadata
#' table that captures library details and FASTQ mapping.
#'
#' The function expects your Pristine Seas Excel workbook with
#' `metadata` and `labwork` sheets. You can pass a FASTQ manifest (data.frame
#' or path to CSV/TSV) with columns like `sample_name`, `fastq_1`, `fastq_2`.
#'
#' @param metadata Path to the Pristine Seas Excel workbook.
#' @param outdir Output directory for the generated TSVs (default: "sra_out").
#' @param fastq_manifest Optional data.frame or CSV/TSV path mapping samples to FASTQ files.
#'
#' @return (invisible) list with paths to the BioSample and SRAmetadata TSVs.
#' @examples
#' \dontrun{
#' pristinesra("expedition_metadata.xlsx",
#'             outdir = "sra_out",
#'             fastq_manifest = "fastq_manifest.tsv")
#' }
#' @importFrom readxl read_excel
#' @importFrom readr write_tsv read_tsv read_csv
#' @importFrom dplyr inner_join mutate transmute rename left_join select filter if_any all_of coalesce cur_data_all
#' @importFrom lubridate parse_date_time
#' @export

pristinesra <- function(metadata, outdir = "sra_out", fastq_manifest = NULL) {
  stopifnot(is.character(metadata), length(metadata) == 1)
  metadata_path <- normalizePath(metadata, mustWork = TRUE)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  ## ---------- helpers (vectorized) ----------
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
  format_date <- function(x) {
    parsed <- suppressWarnings(lubridate::parse_date_time(x, orders = c("ymd","mdy","dmy","ydm")))
    as.character(as.Date(parsed))
  }
  format_latlon <- function(lat, lon) {
    lat_num <- suppressWarnings(as.numeric(lat))
    lon_num <- suppressWarnings(as.numeric(lon))
    out <- rep("", length(lat_num))
    ok  <- is.finite(lat_num) & is.finite(lon_num) &
      lat_num >= -90 & lat_num <= 90 & lon_num >= -180 & lon_num <= 180
    if (any(ok)) {
      lat_str <- formatC(lat_num[ok], format = "f", digits = 8, drop0trailing = TRUE)
      lon_str <- formatC(lon_num[ok], format = "f", digits = 8, drop0trailing = TRUE)
      out[ok] <- paste(lat_str, lon_str)
    }
    out
  }
  to_layout <- function(x) {
    x <- as.character(x); x <- trimws(tolower(x))
    out <- rep(NA_character_, length(x))
    out[x %in% c("paired","paired-end","pair","pe")]  <- "PAIRED"
    out[x %in% c("single","single-end","se")]         <- "SINGLE"
    out
  }
  match_allowed <- function(ans, allowed) {
    m <- pmatch(tolower(ans), tolower(allowed))
    if (is.na(m)) NA_character_ else allowed[m]
  }
  get_choice <- function(prompt_msg, allowed, allow_pass = TRUE) {
    repeat {
      hint <- if (allow_pass) " (type 'pass' to skip, 'h' for help): " else " ('h' for help): "
      ans  <- readline(paste0(prompt_msg, hint))
      ansl <- tolower(trimws(ans))
      if (ansl == "h") {
        cat("\nValid options:\n - ", paste(allowed, collapse = "\n - "), "\n",
            if (allow_pass) "\nType 'pass' to skip.\n\n" else "\n\n", sep = ""); next
      }
      if (allow_pass && ansl == "pass") return("")
      if (ansl == "") { cat("Please enter a value or type 'h' for help",
                            if (allow_pass) ", or 'pass' to skip", ".\n", sep = ""); next }
      fixed <- match_allowed(ans, allowed)
      if (!is.na(fixed)) return(fixed)
      cat("Not recognized. Type 'h' to see valid options",
          if (allow_pass) ", or 'pass' to skip", ".\n", sep = "")
    }
  }
  
  ## ---------- controlled vocabs ----------
  vocab_strategy <- c("WGA","WGS","WXS","RNA-Seq","miRNA-Seq","WCS","CLONE","POOLCLONE","AMPLICON",
                      "CLONEEND","FINISHING","ChIP-Seq","MNase-Seq","DNase-Hypersensitivity",
                      "Bisulfite-Seq","Tn-Seq","EST","FL-cDNA","CTS","MRE-Seq","MeDIP-Seq","MBD-Seq",
                      "Synthetic-Long-Read","ATAC-seq","ChIA-PET","FAIRE-seq","Hi-C","ncRNA-Seq",
                      "RAD-Seq","RIP-Seq","SELEX","ssRNA-seq","Targeted-Capture",
                      "Tethered Chromatin Conformation Capture","DIP-Seq","GBS","Inverse rRNA",
                      "NOMe-Seq","Ribo-seq","VALIDATION","OTHER")
  
  vocab_source <- c("GENOMIC","TRANSCRIPTOMIC","METAGENOMIC","METATRANSCRIPTOMIC","SYNTHETIC",
                    "VIRAL RNA","GENOMIC SINGLE CELL","TRANSCRIPTOMIC SINGLE CELL","OTHER")
  
  vocab_selection <- c("RANDOM","PCR","RANDOM PCR","RT-PCR","HMPR","MF","CF-S","CF-M","CF-H","CF-T","MDA",
                       "MSLL","cDNA","ChIP","MNase","DNAse","Hybrid Selection","Reduced Representation",
                       "Restriction Digest","5-methylcytidine antibody",
                       "MBD2 protein methyl-CpG binding domain","CAGE","RACE","size fractionation",
                       "Padlock probes capture method","other","unspecified","cDNA_oligo_dT",
                       "cDNA_randomPriming","Inverse rRNA","Oligo-dT","PolyA","repeat fractionation")
  
  vocab_platforms <- c("_LS454","ABI_SOLID","BGISEQ","CAPILLARY","COMPLETE_GENOMICS","DNBSEQ","ELEMENT",
                       "GENAPSYS","GENEMIND","HELICOS","ILLUMINA","ION_TORRENT","OXFORD_NANOPORE",
                       "PACBIO_SMRT","TAPESTRI","ULTIMA","VELA_DIAGNOSTICS")
  
  vocab_filetype <- c("bam","srf","sff","fastq","454_native","Helicos_native","SOLiD_native",
                      "PacBio_HDF5","CompleteGenomics_native","OxfordNanopore_native")
  
  model2platform <- c(
    # Illumina
    "Illumina Genome Analyzer"="ILLUMINA","Illumina Genome Analyzer II"="ILLUMINA",
    "Illumina Genome Analyzer IIx"="ILLUMINA","Illumina HiSeq 1000"="ILLUMINA",
    "Illumina HiSeq 1500"="ILLUMINA","Illumina HiSeq 2000"="ILLUMINA",
    "Illumina HiSeq 2500"="ILLUMINA","Illumina HiSeq 3000"="ILLUMINA",
    "Illumina HiSeq 4000"="ILLUMINA","Illumina HiSeq X"="ILLUMINA",
    "Illumina HiSeq X Five"="ILLUMINA","Illumina HiSeq X Ten"="ILLUMINA",
    "Illumina MiSeq"="ILLUMINA","Illumina MiSeq FGx"="ILLUMINA",
    "Illumina MiniSeq"="ILLUMINA","Illumina NovaSeq 6000"="ILLUMINA",
    "Illumina NovaSeq X"="ILLUMINA","Illumina NovaSeq X Plus"="ILLUMINA",
    "Illumina iSeq 100"="ILLUMINA","NextSeq 1000"="ILLUMINA","NextSeq 2000"="ILLUMINA",
    "NextSeq 500"="ILLUMINA","NextSeq 550"="ILLUMINA",
    # Ion Torrent
    "Ion Torrent PGM"="ION_TORRENT","Ion Torrent Proton"="ION_TORRENT",
    "Ion Torrent S5"="ION_TORRENT","Ion Torrent S5 XL"="ION_TORRENT",
    "Ion GeneStudio S5"="ION_TORRENT","Ion GeneStudio S5 Prime"="ION_TORRENT",
    # ONT
    "MinION"="OXFORD_NANOPORE","GridION"="OXFORD_NANOPORE","PromethION"="OXFORD_NANOPORE",
    # PacBio
    "PacBio RS"="PACBIO_SMRT","PacBio RS II"="PACBIO_SMRT","Sequel"="PACBIO_SMRT",
    "Sequel II"="PACBIO_SMRT","Sequel IIe"="PACBIO_SMRT","Onso"="PACBIO_SMRT",
    # 454
    "454 GS"="_LS454","454 GS 20"="_LS454","454 GS FLX"="_LS454","454 GS FLX+"="_LS454",
    "454 GS FLX Titanium"="_LS454","454 GS Junior"="_LS454",
    # ABI / Capillary
    "AB 5500 Genetic Analyzer"="ABI_SOLID","AB 5500xl Genetic Analyzer"="ABI_SOLID",
    "AB SOLiD System"="ABI_SOLID","AB SOLiD 3 Plus System"="ABI_SOLID","AB SOLiD 4 System"="ABI_SOLID",
    "AB SOLiD 4hq System"="ABI_SOLID","AB SOLiD PI System"="ABI_SOLID","AB SOLiD System 2.0"="ABI_SOLID",
    "AB SOLiD System 3.0"="ABI_SOLID",
    "AB 3130 Genetic Analyzer"="CAPILLARY","AB 3130xl Genetic Analyzer"="CAPILLARY",
    "AB 3500 Genetic Analyzer"="CAPILLARY","AB 3500xL Genetic Analyzer"="CAPILLARY",
    "AB 3730 Genetic Analyzer"="CAPILLARY","AB 3730xl Genetic Analyzer"="CAPILLARY",
    # Others
    "BGISEQ-50"="BGISEQ","BGISEQ-500"="BGISEQ","MGISEQ-2000RS"="BGISEQ",
    "DNBSEQ-G400"="DNBSEQ","DNBSEQ-G50"="DNBSEQ","DNBSEQ-T7"="DNBSEQ","DNBSEQ-G400 FAST"="DNBSEQ",
    "Helicos HeliScope"="HELICOS","Complete Genomics"="COMPLETE_GENOMICS",
    "Element AVITI"="ELEMENT","GS111"="GENAPSYS",
    "FASTASeq 300"="GENEMIND","GenoCare 1600"="GENEMIND","GenoLab M"="GENEMIND",
    "Tapestri"="TAPESTRI","UG 100"="ULTIMA","Sentosa SQ301"="VELA_DIAGNOSTICS"
  )
  
  ## ---------- read/clean once ----------
  md <- readxl::read_excel(metadata_path, sheet = "metadata")
  lw <- readxl::read_excel(metadata_path, sheet = "labwork")
  
  for (nm in c("env_local_scale","env_broad_scale","env_medium","Atoll","country","date",
               "depth_m","lat","long","expedition","team_lead"))
    if (!nm %in% names(md)) md[[nm]] <- NA_character_
  for (nm in c("lib_layout","seq_meth","Adaptor Tag","Amplicon Target"))
    if (!nm %in% names(lw)) lw[[nm]] <- NA_character_
  
  md$lat  <- vapply(md$lat,  convert_to_decimal, numeric(1))
  md$long <- vapply(md$long, convert_to_decimal, numeric(1))
  md$date <- format_date(md$date)
  
  ## ---------- interactive SRA type + environment ----------
  cat("\nSelect SRA submission type:\n",
      "1: Environmental metagenomic sequencing (MIMS)\n",
      "2: Metagenome-assembled genome (MIMAG)\n",
      "3: Marker gene sequences (MIMARKS)\n",
      "h: Help describing the options\n", sep = "")
  repeat {
    sra_choice <- readline("Enter 1, 2, 3 or h: ")
    if (tolower(sra_choice) == "h") {
      cat("\nDescriptions:\n",
          "1. MIMS: Environmental/metagenome sequences; scientific name ends with 'metagenome'.\n",
          "2. MIMAG: Metagenome-assembled genomes; no 'metagenome' in organism.\n",
          "3. MIMARKS: Marker gene sequences (e.g., 16S/18S/COI); no 'metagenome' in organism.\n\n", sep = "")
    } else if (sra_choice %in% c("1","2","3")) break
    else cat("Invalid input. Try again.\n")
  }
  sra_type <- switch(sra_choice, "1" = "MIMS", "2" = "MIMAG", "3" = "MIMARKS")
  
  is_ocean  <- tolower(readline("Is this from an ocean environment? (yes/no): "))
  env_input <- if (is_ocean %in% c("no","n")) readline("Enter general environment (e.g., soil, freshwater): ") else "ocean"
  
  ## ---------- BioSample (first) ----------
  metadata_df <- md
  labwork_df  <- lw
  
  merged <- dplyr::inner_join(metadata_df, labwork_df, by = "ps_sample_id") %>%
    dplyr::mutate(
      sample_name        = ps_sample_id,
      sample_title       = "",
      bioproject_accession = "",
      organism           = if (sra_type == "MIMAG") env_input else paste(env_input, "metagenome"),
      isolate            = if (sra_type == "MIMAG") ps_sample_id else "",
      collection_date    = md$date,   # already formatted
      depth              = depth_m,
      env_broad_scale    = env_broad_scale,
      env_local_scale    = env_local_scale,
      env_medium         = env_medium,
      geo_loc_name       = country,
      lat_lon            = format_latlon(md$lat, md$long),
      isolation_source   = paste(env_input, "sample collected in", ifelse(!is.na(Atoll), Atoll, country)),
      metagenome_source  = if (sra_type == "MIMAG") "" else env_input
    )
  
  col_list <- list()
  col_list$MIMS <- c("sample_name","sample_title","bioproject_accession","organism","collection_date","depth",
                     "env_broad_scale","env_local_scale","env_medium","geo_loc_name","lat_lon","alkalinity",
                     "alkalinity_method","alkyl_diethers","altitude","aminopept_act","ammonium","atmospheric_data",
                     "bac_prod","bac_resp","bacteria_carb_prod","biomass","bishomohopanol","bromide","calcium",
                     "carb_nitro_ratio","chem_administration","chloride","chlorophyll","collection_method","conduc",
                     "density","diether_lipids","diss_carb_dioxide","diss_hydrogen","diss_inorg_carb","diss_inorg_nitro",
                     "diss_inorg_phosp","diss_org_carb","diss_org_nitro","diss_oxygen","down_par","elev","fluor",
                     "glucosidase_act","isolation_source","light_intensity","magnesium","mean_frict_vel",
                     "mean_peak_frict_vel","misc_param","n_alkanes","neg_cont_type","nitrate","nitrite","nitro",
                     "omics_observ_id","org_carb","org_matter","org_nitro","organism_count","oxy_stat_samp",
                     "part_org_carb","part_org_nitro","perturbation","petroleum_hydrocarb","ph","phaeopigments",
                     "phosphate","phosplipid_fatt_acid","photon_flux","pos_cont_type","potassium","pressure",
                     "primary_prod","redox_potential","ref_biomaterial","rel_to_oxygen","salinity","samp_collect_device",
                     "samp_mat_process","samp_size","samp_store_dur","samp_store_loc","samp_store_temp",
                     "samp_vol_we_dna_ext","silicate","size_frac","size_frac_low","size_frac_up","sodium",
                     "soluble_react_phosp","source_material_id","sulfate","sulfide","suspend_part_matter","temp",
                     "tidal_stage","tot_depth_water_col","tot_diss_nitro","tot_inorg_nitro","tot_nitro",
                     "tot_part_carb","tot_phosp","turbidity","water_current","description")
  col_list$MIMAG <- c("sample_name","sample_title","bioproject_accession","organism","isolate","collection_date","depth",
                      "env_broad_scale","env_local_scale","env_medium","geo_loc_name","isolation_source","lat_lon",
                      "alkalinity","alkalinity_method","alkyl_diethers","altitude","aminopept_act","ammonium",
                      "atmospheric_data","bac_prod","bac_resp","bacteria_carb_prod","biomass","bishomohopanol",
                      "bromide","calcium","carb_nitro_ratio","chem_administration","chloride","chlorophyll",
                      "collection_method","conduc","density","derived_from","diether_lipids","diss_carb_dioxide",
                      "diss_hydrogen","diss_inorg_carb","diss_inorg_nitro","diss_inorg_phosp","diss_org_carb",
                      "diss_org_nitro","diss_oxygen","down_par","elev","experimental_factor","fluor","glucosidase_act",
                      "light_intensity","magnesium","mean_frict_vel","mean_peak_frict_vel","metagenome_source",
                      "misc_param","n_alkanes","neg_cont_type","nitrate","nitrite","nitro","omics_observ_id",
                      "org_carb","org_matter","org_nitro","organism_count","oxy_stat_samp","part_org_carb",
                      "part_org_nitro","perturbation","petroleum_hydrocarb","ph","phaeopigments","phosphate",
                      "phosplipid_fatt_acid","photon_flux","pos_cont_type","potassium","pressure","primary_prod",
                      "redox_potential","ref_biomaterial","rel_to_oxygen","salinity","samp_collect_device",
                      "samp_mat_process","samp_size","samp_store_dur","samp_store_loc","samp_store_temp",
                      "samp_vol_we_dna_ext","silicate","size_frac","size_frac_low","size_frac_up","sodium",
                      "soluble_react_phosp","source_material_id","sulfate","sulfide","suspend_part_matter","temp",
                      "tidal_stage","tot_depth_water_col","tot_diss_nitro","tot_inorg_nitro","tot_nitro",
                      "tot_part_carb","tot_phosp","turbidity","water_current","description")
  col_list$MIMARKS <- setdiff(col_list$MIMS, "isolate")
  
  miss_cols <- setdiff(col_list[[sra_type]], names(merged))
  for (nm in miss_cols) merged[[nm]] <- ""
  bio_df  <- merged[, col_list[[sra_type]], drop = FALSE]
  bio_path <- file.path(outdir, paste0("BioSample_", tolower(sra_type), ".tsv"))
  readr::write_tsv(bio_df, bio_path, na = "")
  message("BioSample file written to: ", normalizePath(bio_path, winslash = "/"))
  
  ## ---------- SRAmetadata (second) ----------
  df <- dplyr::inner_join(md, lw, by = "ps_sample_id") %>%
    dplyr::mutate(
      sample_name        = ps_sample_id,
      library_ID         = if ("Adaptor Tag" %in% names(.)) paste0(ps_sample_id, "_", `Adaptor Tag`) else ps_sample_id,
      title              = paste0("eDNA ", ps_sample_id),
      library_layout     = to_layout(dplyr::coalesce(lib_layout, NA_character_)),
      design_description = dplyr::case_when(
        !is.na(`Amplicon Target`) & `Amplicon Target` != "" ~ paste("Marker:", `Amplicon Target`),
        TRUE ~ "eDNA sequencing"
      ),
      instrument_model   = if ("seq_meth" %in% names(.)) dplyr::coalesce(seq_meth, NA_character_) else NA_character_,
      platform           = NA_character_
    )
  
  for (nm in c("library_source","library_selection","library_strategy",
               "filetype","filename","filename2","filename3","filename4","assembly","fasta_file")) {
    if (!nm %in% names(df)) df[[nm]] <- NA_character_
  }
  
  if (!is.null(fastq_manifest)) {
    man <- if (is.data.frame(fastq_manifest)) {
      fastq_manifest
    } else if (is.character(fastq_manifest) && file.exists(fastq_manifest)) {
      if (grepl("\\.csv$", fastq_manifest, TRUE))
        readr::read_csv(fastq_manifest, show_col_types = FALSE)
      else
        readr::read_tsv(fastq_manifest, show_col_types = FALSE)
    } else stop("`fastq_manifest` must be a data.frame or a readable CSV/TSV path.")
    names(man) <- tolower(names(man))
    if (!"sample_name" %in% names(man)) {
      alt <- intersect(c("sampleid","sample_id","sample"), names(man))
      if (length(alt)) man <- dplyr::rename(man, sample_name = .data[[alt[1]]])
      else stop("FASTQ manifest needs a `sample_name` column.")
    }
    pick <- function(df2, cand) { nm <- cand[cand %in% names(df2)][1]; if (length(nm)) df2[[nm]] else NULL }
    man <- dplyr::transmute(
      man,
      sample_name,
      filetype  = pick(dplyr::cur_data_all(), c("filetype")),
      filename  = pick(dplyr::cur_data_all(), c("fastq_1","fastq1","r1","read1","file1","forward")),
      filename2 = pick(dplyr::cur_data_all(), c("fastq_2","fastq2","r2","read2","file2","reverse")),
      filename3 = pick(dplyr::cur_data_all(), c("fastq_3","file3")),
      filename4 = pick(dplyr::cur_data_all(), c("fastq_4","file4"))
    )
    df <- dplyr::left_join(df, man, by = "sample_name", suffix = c("", ".man"))
    for (nm in c("filetype","filename","filename2","filename3","filename4"))
      df[[nm]] <- dplyr::coalesce(df[[paste0(nm, ".man")]], df[[nm]])
    df <- dplyr::select(df, -dplyr::ends_with(".man"))
  }
  
  idx <- is.na(df$platform) & !is.na(df$instrument_model)
  if (any(idx)) df$platform[idx] <- unname(model2platform[df$instrument_model[idx]])
  
  if (all(is.na(df$instrument_model)) && all(is.na(df$platform))) {
    im <- get_choice("Enter instrument_model", allowed = c(names(model2platform)), allow_pass = TRUE)
    if (im != "") {
      df$instrument_model[is.na(df$instrument_model)] <- im
      df$platform[is.na(df$platform)] <- unname(model2platform[im])
    } else {
      pl <- get_choice("Enter platform", allowed = vocab_platforms, allow_pass = TRUE)
      if (pl != "") df$platform[is.na(df$platform)] <- pl
    }
  }
  
  if (all(is.na(df$library_source))) {
    ls <- get_choice("Choose library_source", allowed = vocab_source, allow_pass = TRUE)
    if (ls != "") df$library_source[is.na(df$library_source) | df$library_source == ""] <- ls
  }
  if (all(is.na(df$library_selection))) {
    sel <- get_choice("Choose library_selection", allowed = vocab_selection, allow_pass = TRUE)
    if (sel != "") df$library_selection[is.na(df$library_selection) | df$library_selection == ""] <- sel
  }
  if (all(is.na(df$library_strategy))) {
    st <- get_choice("Choose library_strategy", allowed = vocab_strategy, allow_pass = TRUE)
    if (st != "") df$library_strategy[is.na(df$library_strategy) | df$library_strategy == ""] <- st
  }
  if (all(is.na(df$filetype))) {
    ft <- get_choice("Choose filetype", allowed = vocab_filetype, allow_pass = TRUE)
    if (ft != "") df$filetype[is.na(df$filetype) | df$filetype == ""] <- ft
  }
  
  if (all(is.na(df$library_layout))) {
    ll <- get_choice("Set default library_layout (PAIRED/SINGLE)", allowed = c("PAIRED","SINGLE"), allow_pass = FALSE)
    df$library_layout <- ll
  } else if (any(is.na(df$library_layout))) {
    ans <- get_choice("Set default library_layout for missing rows (PAIRED/SINGLE)", allowed = c("PAIRED","SINGLE"), allow_pass = FALSE)
    df$library_layout[is.na(df$library_layout)] <- ans
  }
  
  req_sra <- c("sample_name","library_ID","title","library_layout","design_description")
  miss_mat <- sapply(req_sra, function(cn) is.na(df[[cn]]) | df[[cn]] == "")
  if (is.matrix(miss_mat) && any(miss_mat)) {
    miss_tbl <- dplyr::transmute(df, sample_name, !!!setNames(as.data.frame(miss_mat), req_sra)) %>%
      dplyr::filter(dplyr::if_any(dplyr::all_of(req_sra), ~ .x))
    if (nrow(miss_tbl)) {
      miss_path <- file.path(outdir, "SRA_required_missing.tsv")
      readr::write_tsv(miss_tbl, miss_path, na = "")
      message("Missing required SRA fields written to ", normalizePath(miss_path, winslash = "/"))
    }
  }
  
  out_cols <- c("sample_name","library_ID","title","library_strategy","library_source","library_selection",
                "library_layout","platform","instrument_model","design_description",
                "filetype","filename","filename2","filename3","filename4","assembly","fasta_file")
  for (nm in setdiff(out_cols, names(df))) df[[nm]] <- NA_character_
  sra_md   <- df[, out_cols, drop = FALSE]
  sra_path <- file.path(outdir, "SRAmetadata.tsv")
  readr::write_tsv(sra_md, sra_path, na = "")
  message("Wrote: ", normalizePath(sra_path, winslash = "/"))
  
  invisible(list(biosample = bio_path, sra = sra_path))
}
