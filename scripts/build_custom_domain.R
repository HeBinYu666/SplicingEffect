#!/usr/bin/env Rscript

library(AnnotationHub)
library(ensembldb)
library(GenomicRanges)
library(dplyr)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/build_custom_domain.R --pfam pfam_results.txt --out domain.csv --ensembl 110\n\n",
    "Options:\n",
    "  --pfam      PfamScan result table\n",
    "  --out       Output domain-coordinate CSV used by IsoImpact with -d/--domain\n",
    "  --ensembl   Ensembl release version number for the AnnotationHub EnsDb query. Default: 110\n",
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

# 1. 读取用户指定的输入和输出路径
pfam_txt <- get_arg("--pfam")
output_csv <- get_arg("--out")
ensembl_release <- get_arg("--ensembl", "110")
species <- get_arg("--species", "Homo sapiens")
gtf_file <- get_arg("--gtf")

if (is.null(pfam_txt) || is.null(output_csv)) {
  usage()
  stop("必须提供 --pfam 和 --out 参数。", call. = FALSE)
}

if (!file.exists(pfam_txt)) {
  stop(paste("报错：找不到 Pfam 文件：", pfam_txt), call. = FALSE)
}

if (!is.null(gtf_file)) {
  print(paste("已指定匹配的 GTF 文件，用于后续 IsoImpact 分析：", gtf_file))
}

# 2. 按照官方 AnnotationHub 流程获取指定 Ensembl 版本的 EnsDb 数据库
print(paste("正在连接 AnnotationHub 获取 Ensembl v", ensembl_release, " 数据库...", sep = ""))
ah <- AnnotationHub()
query_result <- query(ah, c("EnsDb", species, ensembl_release))

if (length(query_result) == 0) {
  stop(
    paste(
      "报错：AnnotationHub 中没有找到",
      species,
      "Ensembl v",
      ensembl_release,
      "对应的 EnsDb 数据库。",
      sep = " "
    ),
    call. = FALSE
  )
}

edb <- query_result[[1]]

# 3. 严谨读取 PfamScan 结果（跳过 # 注释行，自动分列）
print("正在读取 Pfam 结果...")
pfam_data <- read.table(pfam_txt, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)

if (nrow(pfam_data) == 0) {
  stop("报错：Pfam 文件是空的或者没有正确读取，请检查输入文件！", call. = FALSE)
}

if (ncol(pfam_data) < 7) {
  stop("报错：Pfam 文件列数少于 7 列，请确认是否为标准 PfamScan 输出。", call. = FALSE)
}

print(paste("成功读取了", nrow(pfam_data), "条 Pfam Domain 记录，开始映射..."))

results <- list()
valid_count <- 0

# 4. 批量执行映射
for (i in 1:nrow(pfam_data)) {
  row <- pfam_data[i, ]

  # pfam_scan.pl 的标准输出列：V1 是 ID，V2 是起始，V3 是终止，V6 是 Pfam ID，V7 是 Domain 名字
  prot_id_raw <- as.character(row$V1)
  # 必须切掉小数点版本号，否则 proteinToGenome 不认
  protein_id <- strsplit(prot_id_raw, "\\.")[[1]][1]

  pfam_start <- as.numeric(row$V2)
  pfam_end <- as.numeric(row$V3)
  domain_name <- as.character(row$V7)

  # 提取 V6 列的 Pfam ID，例如 PF00022.20
  pfam_id_raw <- as.character(row$V6)
  # 把小数点及后面的版本号去掉，只保留 PFxxxx
  domain_id <- strsplit(pfam_id_raw, "\\.")[[1]][1]

  tryCatch({
    protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id)
    genomic_coords <- proteinToGenome(protein_range, db = edb)

    if (length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
      g_start <- min(start(genomic_coords[[1]]))
      g_end <- max(end(genomic_coords[[1]]))

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
    # 遇到无法映射的孤儿序列静默跳过
  })

  if (i %% 1000 == 0) {
    print(paste("已处理", i, "条记录... 当前成功映射:", valid_count, "条"))
  }
}

# 5. 输出 CSV
if (valid_count > 0) {
  final_df <- bind_rows(results)
  final_df <- distinct(final_df, Protein_ID, Domain_ID, Genomic_Start, Genomic_End, .keep_all = TRUE)
  write.csv(final_df, output_csv, row.names = FALSE)
  print(paste("全量绝对坐标数据库构建完成！共成功映射了", nrow(final_df), "个 Domain。"))
} else {
  print("警告：没有成功映射任何 Domain，请检查 ID 是否匹配。")
}
