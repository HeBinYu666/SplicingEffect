#!/usr/bin/env Rscript

# Build an IsoImpact-compatible domain-coordinate CSV using the official
# AnnotationHub/EnsDb protein-to-genome mapping workflow.
# Required input:
#   1. A standard PfamScan result table. The script expects protein ID,
#      domain start, domain end, Pfam ID, and domain name in the standard columns.
#   2. Either an Ensembl release number with a matching EnsDb record in
#      AnnotationHub, or a local EnsDb sqlite file.

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
  library(GenomicRanges)
  library(dplyr)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_custom_domain.R --pfam pfam_results.txt --out domain.csv --ensembl 110\n\n",
    "Options:\n",
    "  --pfam      PfamScan result table for the corresponding protein sequences\n",
    "  --out       Output CSV used by IsoImpact with -d/--domain\n",
    "  --ensembl   Ensembl release number for AnnotationHub EnsDb lookup. Default: 110\n",
    "  --species   Species name used in AnnotationHub query. Default: Homo sapiens\n",
    "  --ensdb-sqlite  Optional local EnsDb sqlite file with protein annotations. If supplied, AnnotationHub is not used\n",
    "  --gtf       Optional matching GTF path. Accepted for workflow consistency; not used for mapping\n",
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

PFAM_TXT <- get_arg("--pfam")
OUTPUT_CSV <- get_arg("--out", "custom_domain.csv")
ENSEMBL_RELEASE <- get_arg("--ensembl", "110")
SPECIES <- get_arg("--species", "Homo sapiens")
USER_GTF <- get_arg("--gtf")
ENSDB_SQLITE <- get_arg("--ensdb-sqlite")

if (is.null(PFAM_TXT)) {
  usage()
  stop("--pfam is required.", call. = FALSE)
}

if (!is.null(USER_GTF)) {
  message("[IsoImpact] Matching GTF supplied for the downstream IsoImpact run: ", USER_GTF)
}

if (!is.null(ENSDB_SQLITE)) {
  message("[IsoImpact] Loading local EnsDb sqlite file: ", ENSDB_SQLITE)
  if (!file.exists(ENSDB_SQLITE)) stop("EnsDb sqlite file not found: ", ENSDB_SQLITE, call. = FALSE)
  edb <- EnsDb(ENSDB_SQLITE)
} else {
  message("[IsoImpact] Loading EnsDb from AnnotationHub...")
  message("[IsoImpact] Species: ", SPECIES)
  message("[IsoImpact] Ensembl release: ", ENSEMBL_RELEASE)

  ah <- AnnotationHub()
  query_result <- query(ah, c("EnsDb", SPECIES, ENSEMBL_RELEASE))
  if (length(query_result) == 0) {
    stop(
      "No matching EnsDb record was found in AnnotationHub for species '",
      SPECIES, "' and Ensembl release '", ENSEMBL_RELEASE, "'.",
      call. = FALSE
    )
  }
  edb <- query_result[[1]]
  message("[IsoImpact] EnsDb loaded: ", names(query_result)[1])
}

message("[IsoImpact] Reading PfamScan results...")
if (!file.exists(PFAM_TXT)) stop("PfamScan result file not found: ", PFAM_TXT, call. = FALSE)

pfam_data <- read.table(PFAM_TXT, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
if (nrow(pfam_data) == 0) stop("PfamScan result file is empty or could not be parsed.", call. = FALSE)
if (ncol(pfam_data) < 7) {
  stop(
    "PfamScan result file has fewer than seven columns. Expected standard PfamScan output.",
    call. = FALSE
  )
}

message("[IsoImpact] PfamScan records detected: ", nrow(pfam_data))

results <- list()
valid_count <- 0L

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
      genomic_coords <- proteinToGenome(protein_range, db = edb)
    }, error = function(e) {})

    if (is.null(genomic_coords) || length(genomic_coords) == 0) {
      tryCatch({
        protein_range <- IRanges(start = pfam_start, end = pfam_end, names = prot_id_raw)
        genomic_coords <- proteinToGenome(protein_range, db = edb)
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

    valid_count <- valid_count + 1L
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
    message("[IsoImpact] Processed ", i, " PfamScan records; mapped domains: ", valid_count)
  }
}

if (valid_count > 0) {
  final_df <- bind_rows(results)
  final_df <- distinct(final_df, Protein_ID, Domain_ID, Genomic_Start, Genomic_End, .keep_all = TRUE)

  full_csv_name <- sub("\\.csv$", "_full.csv", OUTPUT_CSV)
  write.csv(final_df, full_csv_name, row.names = FALSE)

  final_df_compat <- final_df[, c("Protein_ID", "Domain_ID", "Domain_Name", "Genomic_Start", "Genomic_End")]
  write.csv(final_df_compat, OUTPUT_CSV, row.names = FALSE)

  message("[IsoImpact] Domain CSV written to: ", OUTPUT_CSV)
  message("[IsoImpact] Full mapping table written to: ", full_csv_name)
  message("[IsoImpact] Mapped non-redundant domains: ", nrow(final_df_compat))
} else {
  warning("No PfamScan domains could be mapped to genomic coordinates.")
}
