#!/usr/bin/env Rscript

# Build an IsoImpact-compatible domain-coordinate CSV with official
# ensembldb coordinate-mapping functions.
# Required inputs:
#   1. A CDS-aware GTF file whose protein_id values match the protein IDs used
#      in the Pfam/PfamScan results.
#   2. A standard Pfam/PfamScan result table. The script expects protein ID,
#      domain start, domain end, Pfam ID, and domain name in the standard columns.

suppressPackageStartupMessages({
  library(ensembldb)
  library(GenomicRanges)
  library(dplyr)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_custom_domain.R --gtf annotation.gtf --pfam pfam_results.txt --out domain.csv\n\n",
    "Options:\n",
    "  --gtf      CDS-aware GTF file for the target annotation version\n",
    "  --pfam     Pfam/PfamScan result table for the corresponding protein sequences\n",
    "  --out      Output CSV used by IsoImpact with -d/--domain\n",
    "  --sqlite   Optional temporary EnsDb sqlite file path\n",
    "  --help     Show this message\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || length(args) == 0) {
  usage()
  quit(status = 0)
}

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag), call. = FALSE)
  args[[idx + 1]]
}

USER_GTF <- get_arg("--gtf")
PFAM_TXT <- get_arg("--pfam")
OUTPUT_CSV <- get_arg("--out", "custom_domain.csv")
CUSTOM_SQLITE <- get_arg("--sqlite", sub("\\.csv$", ".sqlite", OUTPUT_CSV))

if (is.null(USER_GTF) || is.null(PFAM_TXT)) {
  usage()
  stop("--gtf and --pfam are required.", call. = FALSE)
}

message("[IsoImpact] Checking custom GTF...")
if (!file.exists(USER_GTF)) stop("GTF file not found: ", USER_GTF, call. = FALSE)

open_gtf <- if (grepl("\\.gz$", USER_GTF, ignore.case = TRUE)) gzfile else file
con <- open_gtf(USER_GTF, open = "rt")
on.exit(close(con), add = TRUE)

cds_count <- 0L
repeat {
  lines <- readLines(con, n = 100000, warn = FALSE)
  if (!length(lines)) break
  cds_count <- cds_count + sum(grepl("\tCDS\t", lines, fixed = TRUE))
}

if (is.na(cds_count) || cds_count == 0) {
  stop(
    "The GTF file does not contain CDS records. proteinToGenome requires CDS ",
    "coordinates. Add CDS annotation first, for example with TransDecoder or gffread.",
    call. = FALSE
  )
}
message("[IsoImpact] CDS records detected: ", cds_count)

message("[IsoImpact] Building temporary EnsDb database...")
if (file.exists(CUSTOM_SQLITE)) file.remove(CUSTOM_SQLITE)

suppressMessages({
  ensDbFromGtf(
    USER_GTF,
    outfile = CUSTOM_SQLITE,
    organism = "Custom_species",
    genomeVersion = "Custom_v1",
    version = 1
  )
})
custom_db <- EnsDb(CUSTOM_SQLITE)

message("[IsoImpact] Reading Pfam results...")
if (!file.exists(PFAM_TXT)) stop("Pfam result file not found: ", PFAM_TXT, call. = FALSE)

pfam_data <- read.table(PFAM_TXT, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
if (nrow(pfam_data) == 0) stop("Pfam result file is empty or could not be parsed.", call. = FALSE)
if (ncol(pfam_data) < 7) {
  stop(
    "Pfam result file has fewer than seven columns. Expected standard Pfam/PfamScan output.",
    call. = FALSE
  )
}

message("[IsoImpact] Pfam records detected: ", nrow(pfam_data))

results <- list()
valid_count <- 0

for (i in seq_len(nrow(pfam_data))) {
  row <- pfam_data[i, ]

  prot_id_raw <- as.character(row$V1)
  protein_id_no_version <- sub("\\.[0-9]+$", "", prot_id_raw)

  pfam_start <- suppressWarnings(as.numeric(row$V2))
  pfam_end <- suppressWarnings(as.numeric(row$V3))
  if (is.na(pfam_start) || is.na(pfam_end)) next

  domain_id <- sub("\\.[0-9]+$", "", as.character(row$V6))
  domain_name <- as.character(row$V7)
  if (is.na(domain_name) || domain_name == "") domain_name <- domain_id

  genomic_coords <- NULL

  suppressWarnings(suppressMessages({
    tryCatch({
      protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id_no_version)
      genomic_coords <- proteinToGenome(protein_range, db = custom_db)
    }, error = function(e) {})

    if (is.null(genomic_coords) || length(genomic_coords) == 0) {
      tryCatch({
        protein_range <- IRanges(start = pfam_start, end = pfam_end, names = prot_id_raw)
        genomic_coords <- proteinToGenome(protein_range, db = custom_db)
      }, error = function(e) {})
    }
  }))

  if (!is.null(genomic_coords) && length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
    gr <- genomic_coords[[1]]

    chrom <- as.character(seqnames(gr)[1])
    strnd <- as.character(strand(gr)[1])
    g_start <- min(start(gr))
    g_end <- max(end(gr))
    blocks <- paste(start(gr), end(gr), sep = "-", collapse = ";")

    valid_count <- valid_count + 1
    results[[valid_count]] <- data.frame(
      Protein_ID = protein_id_no_version,
      Domain_ID = domain_id,
      Domain_Name = domain_name,
      Chr = chrom,
      Strand = strnd,
      Genomic_Start = g_start,
      Genomic_End = g_end,
      Exon_Blocks = blocks,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 1000 == 0) {
    message("[IsoImpact] Processed ", i, " Pfam records; mapped domains: ", valid_count)
  }
}

if (valid_count > 0) {
  final_df <- bind_rows(results)
  final_df <- distinct(final_df, Protein_ID, Domain_ID, Genomic_Start, Genomic_End, .keep_all = TRUE)

  full_csv_name <- sub("\\.csv$", "_full.csv", OUTPUT_CSV)
  write.csv(final_df, full_csv_name, row.names = FALSE)

  final_df_compat <- final_df[, c("Protein_ID", "Domain_ID", "Domain_Name", "Genomic_Start", "Genomic_End")]
  write.csv(final_df_compat, OUTPUT_CSV, row.names = FALSE)

  message("[IsoImpact] Custom domain CSV written to: ", OUTPUT_CSV)
  message("[IsoImpact] Full mapping table written to: ", full_csv_name)
  message("[IsoImpact] Mapped non-redundant domains: ", nrow(final_df_compat))
} else {
  warning("No Pfam domains could be mapped to genomic coordinates.")
}
