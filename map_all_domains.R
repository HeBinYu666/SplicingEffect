
library(AnnotationHub)
library(ensembldb)
library(GenomicRanges)
library(dplyr)

# 1. 严格使用你指定的路径
pfam_txt <- " " 
output_csv <- " "

# 2. 完全复刻老师的 AnnotationHub 获取流程
print("正在连接 AnnotationHub 获取 Ensembl v110 数据库...")
ah <- AnnotationHub()
query_result <- query(ah, c("EnsDb", "Homo sapiens", "110"))
edb_v110 <- query_result[[1]]

# 3. 严谨读取 Pfam 结果 (跳过 # 注释行，自动分列)
print("正在读取 Pfam 结果...")
pfam_data <- read.table(pfam_txt, stringsAsFactors=FALSE, comment.char="#", fill=TRUE)

if(nrow(pfam_data) == 0) {
  stop("报错：Pfam 文件是空的或者没有正确读取，请检查 all_pfam_results.txt！")
}
print(paste("成功读取了", nrow(pfam_data), "条 Pfam Domain 记录，开始映射..."))

results <- list()
valid_count <- 0

# 4. 批量执行映射
for(i in 1:nrow(pfam_data)) {
  row <- pfam_data[i, ]
  
  # pfam_scan.pl 的标准输出列：V1是ID, V2是起始, V3是终止, V7是Domain名字
  prot_id_raw <- as.character(row$V1)
  # 必须切掉小数点版本号，否则 proteinToGenome 不认 (如 ENSP00000364361.1 -> ENSP00000364361)
  protein_id <- strsplit(prot_id_raw, "\\.")[[1]][1] 
  
  pfam_start <- as.numeric(row$V2)
  pfam_end <- as.numeric(row$V3)
  domain_name <- as.character(row$V7)
  
  tryCatch({
    # 完全复刻你的代码：names = protein_id
    protein_range <- IRanges(start = pfam_start, end = pfam_end, names = protein_id)
    genomic_coords <- proteinToGenome(protein_range, db = edb_v110)
    
    if(length(genomic_coords) > 0 && length(genomic_coords[[1]]) > 0) {
      g_start <- min(start(genomic_coords[[1]]))
      g_end <- max(end(genomic_coords[[1]]))
      
      valid_count <- valid_count + 1
      results[[valid_count]] <- data.frame(
        Protein_ID = protein_id,
        Domain_Name = domain_name,
        Genomic_Start = g_start,
        Genomic_End = g_end,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
      # 遇到无法映射的孤儿序列静默跳过，保证十万条不中断
  })
  
  if(i %% 1000 == 0) print(paste("已处理", i, "条记录... 当前成功映射:", valid_count, "条"))
}

# 5. 输出 CSV
if(valid_count > 0) {
  final_df <- bind_rows(results)
  write.csv(final_df, output_csv, row.names = FALSE)
  print(paste("全量绝对坐标数据库构建完成！共成功映射了", valid_count, "个 Domain。"))
} else {
  print("警告：没有成功映射任何 Domain，请检查 ID 是否匹配。")
}
