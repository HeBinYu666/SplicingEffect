#!/usr/bin/env Rscript

# Reproduce an IsoImpact-compatible reference domain-coordinate file from
# Pfam/PfamScan protein-domain results and an Ensembl EnsDb annotation.
# This script documents how the built-in human/mouse files were generated.

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
  library(GenomicRanges)
  library(dplyr)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_reference_domain_database.R --species \"Homo sapiens\" --release 110 --pfam all_pfam_results.txt --out human_domain.csv\n\n",
    "Options:\n",
    "  --species   Ensembl species name used by AnnotationHub, e.g. \"Homo sapiens\" or \"Mus musculus\"\n",
    "  --release   Ensembl release number, e.g. 110\n",
    "  --pfam      Pfam/PfamScan result table\n",
    "  --out       Output IsoImpact-compatible domain CSV\n",
    "  --help      Show this message\n",
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

SPECIES <- get_arg("--species")
RELEASE <- get_arg("--release")
PFAM_TXT <- get_arg("--pfam")
OUTPUT_CSV <- get_arg("--out", "human_domain.csv")

if (is.null(SPECIES) || is.null(RELEASE) || is.null(PFAM_TXT)) {
  usage()
  stop("--species, --release, and --pfam are required.", call. = FALSE)
}

if (!file.exists(PFAM_TXT)) stop("Pfam result file not found: ", PFAM_TXT, call. = FALSE)

message("[IsoImpact] Querying AnnotationHub for EnsDb: ", SPECIES, " release ", RELEASE)
ah <- AnnotationHub()
query_result <- query(ah, c("EnsDb", SPECIES, as.character(RELEASE)))
if (length(query_result) == 0) {
  stop("No matching EnsDb record found in AnnotationHub.", call. = FALSE)
}
edb <- query_result[[1]]

message("[IsoImpact] Reading Pfam results...")
pfam_data <- read.table(PFAM_TXT, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
if (nrow(pfam_data) == 0) stop("Pfam result file is empty or could not be parsed.", call. = FALSE)
if (ncol(pfam_data) < 7) {
  stop("Pfam result file has fewer than seven columns. Expected standard Pfam/PfamScan output.", call. = FALSE)
}

results <- list()
valid_count <- 0

for (i in seq_len(nrow(pfam_data))) {
  row <- pfam_data[i, ]

  prot_id_raw <- as.character(row$V1)
  protein_id <- sub("\\.[0-9]+$", "", prot_id_raw)
  pfam_start <- suppressWarnings(as.numeric(row$V2))
  pfam_end <- suppressWarnings(as.numeric(row$V3))
  if (is.na(pfam_start) || is.na(pfam_end)) next

  domain_id <- sub("\\.[0-9]+$", "", as.character(row$V6))
  domain_name <- as.character(row$V7)
  if (is.na(domain_name) || domain_name == "") domain_name <- domain_id

  genomic_coords <- NULL
  suppressWarnings(suppressMessages({
    tryCatch({
      protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id)
      genomic_coords <- proteinToGenome(protein_range, db = edb)
    }, error = function(e) {})
  }))

  if (!is.null(genomic_coords) && length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
    gr <- genomic_coords[[1]]

    valid_count <- valid_count + 1
    results[[valid_count]] <- data.frame(
      Protein_ID = protein_id,
      Domain_ID = domain_id,
      Domain_Name = domain_name,
      Genomic_Start = min(start(gr)),
      Genomic_End = max(end(gr)),
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
  write.csv(final_df, OUTPUT_CSV, row.names = FALSE)

  message("[IsoImpact] Reference domain CSV written to: ", OUTPUT_CSV)
  message("[IsoImpact] Mapped non-redundant domains: ", nrow(final_df))
} else {
  warning("No Pfam domains could be mapped to genomic coordinates.")
}
