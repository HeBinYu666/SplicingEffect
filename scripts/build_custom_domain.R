#!/usr/bin/env Rscript

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
    "  --pfam      PfamScan result table\n",
    "  --out       Output domain-coordinate CSV used by IsoImpact with -d/--domain\n",
    "  --ensembl   Ensembl release number for the AnnotationHub EnsDb query. Default: 110\n",
    "  --species   Species name for the AnnotationHub EnsDb query. Default: Homo sapiens\n",
    "  --gtf       Optional matching GTF path for record keeping; the mapping uses AnnotationHub EnsDb\n",
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

pfam_txt <- get_arg("--pfam")
output_csv <- get_arg("--out", "custom_domain.csv")
ensembl_release <- get_arg("--ensembl", "110")
species <- get_arg("--species", "Homo sapiens")
gtf_file <- get_arg("--gtf")

if (is.null(pfam_txt)) {
  usage()
  stop("--pfam is required.", call. = FALSE)
}

if (!file.exists(pfam_txt)) stop("PfamScan result file not found: ", pfam_txt, call. = FALSE)
if (!is.null(gtf_file)) {
  message("[IsoImpact] Matching GTF supplied for the downstream IsoImpact run: ", gtf_file)
}

print(paste("Connecting to AnnotationHub and retrieving Ensembl release", ensembl_release, "EnsDb..."))
ah <- AnnotationHub()
query_result <- query(ah, c("EnsDb", species, ensembl_release))

if (length(query_result) == 0) {
  stop(
    "No EnsDb record was found in AnnotationHub for species '",
    species, "' and Ensembl release '", ensembl_release, "'.",
    call. = FALSE
  )
}

edb <- query_result[[1]]
print(paste("Successfully loaded EnsDb:", names(query_result)[1]))

print("Reading PfamScan results...")
pfam_data <- read.table(pfam_txt, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)

if (nrow(pfam_data) == 0) {
  stop("PfamScan result file is empty or could not be parsed.", call. = FALSE)
}
if (ncol(pfam_data) < 7) {
  stop("PfamScan result file has fewer than seven columns.", call. = FALSE)
}

print(paste("Successfully read", nrow(pfam_data), "PfamScan domain records. Starting coordinate mapping..."))

results <- list()
valid_count <- 0

for (i in seq_len(nrow(pfam_data))) {
  row <- pfam_data[i, ]

  prot_id_raw <- as.character(row$V1)
  protein_id <- strsplit(prot_id_raw, "\\.")[[1]][1]

  pfam_start <- suppressWarnings(as.numeric(row$V2))
  pfam_end <- suppressWarnings(as.numeric(row$V3))
  if (is.na(pfam_start) || is.na(pfam_end)) next

  pfam_id_raw <- as.character(row$V6)
  domain_id <- strsplit(pfam_id_raw, "\\.")[[1]][1]
  domain_name <- as.character(row$V7)
  if (is.na(domain_name) || domain_name == "") domain_name <- domain_id

  tryCatch({
    protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id)
    genomic_coords <- proteinToGenome(protein_range, db = edb)

    if (length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
      gr <- genomic_coords[[1]]
      g_start <- min(start(gr))
      g_end <- max(end(gr))

      valid_count <- valid_count + 1
      results[[valid_count]] <- data.frame(
        Protein_ID = protein_id,
        Domain_ID = domain_id,
        Domain_Name = domain_name,
        Genomic_Start = g_start,
        Genomic_End = g_end,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    # Skip proteins that cannot be mapped by the selected EnsDb release.
  })

  if (i %% 1000 == 0) {
    print(paste("Processed", i, "records; successfully mapped", valid_count, "domains."))
  }
}

if (valid_count > 0) {
  final_df <- bind_rows(results)
  final_df <- distinct(final_df, Protein_ID, Domain_ID, Genomic_Start, Genomic_End, .keep_all = TRUE)

  write.csv(final_df, output_csv, row.names = FALSE)
  print(paste("Domain-coordinate CSV written to:", output_csv))
  print(paste("Mapped non-redundant domains:", nrow(final_df)))
} else {
  warning("No PfamScan domains could be mapped to genomic coordinates.")
}
