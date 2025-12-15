############# SCENT Analysis of zebrafish embryo ISSAAC-seq dataset
# devtools::install_github("immunogenomics/SCENT")
library(SCENT)
library(Signac)
library(Seurat)
library(dplyr)
library(TFBSTools)

########### Create gene-peaks pair datasets(peaks.info)
#### Read genes annotation
library(rtracklayer)
library(GenomicRanges)

# 读取GTF文件
gtf <- import("~/reference/chr_Danio_rerio.GRCz11.112.chr_filtered.gtf")

# 筛选gene条目
genes <- gtf[gtf$type == "gene"]

# 创建包含必要信息的GRanges对象
gene_ranges <- GRanges(
  seqnames = seqnames(genes),
  ranges = ranges(genes),
  strand = strand(genes),
  gene_id = genes$gene_id,
  gene_name = ifelse(is.na(genes$gene_name), genes$gene_id, genes$gene_name)
)

### Read Topic peaks
## TBX
peaks_topic8 <- read.table(file = "~/Multiomics-ZBF/250106_issaac/script/scenicplus/outs/region_sets/Topics_otsu_peaksid/Topic8.bed")
peaks_topic47 <- read.table(file = "~/Multiomics-ZBF/250106_issaac/script/scenicplus/outs/region_sets/Topics_otsu_peaksid/Topic47.bed")
peaks_topic49 <- read.table(file = "~/Multiomics-ZBF/250106_issaac/script/scenicplus/outs/region_sets/Topics_otsu_peaksid/Topic49.bed")

peaks_topic8$peaks <- paste(peaks_topic8$V1 , peaks_topic8$V2+1 , peaks_topic8$V3 , sep = "-")
peaks_topic47$peaks <- paste(peaks_topic47$V1 , peaks_topic47$V2+1 , peaks_topic47$V3 , sep = "-")
peaks_topic49$peaks <- paste(peaks_topic49$V1 , peaks_topic49$V2+1 , peaks_topic49$V3 , sep = "-")

# Create the gene-peak parir of Topic8.
peaks_gr <- StringToGRanges(peaks_topic49$peaks)
genes_gr <- gene_ranges
peaks_centers <- resize(peaks_gr, width = 1, fix = "center")

hits <- findOverlaps(
  query = genes_gr,
  subject = peaks_centers,
  maxgap = 500000,
  ignore.strand = TRUE
)
peaks_info <- data.frame(Gene_name = genes_gr$gene_name[hits@from] , Peak = GRangesToString(peaks_gr[hits@to]))

peaks_info$Gene_name[peaks_info$Gene_name == c("Metazoa_SRP")] <- "Metazoa-SRP"
peaks_info$Gene_name[peaks_info$Gene_name == c("RNase_MRP")] <- "RNase-MRP"
peaks_info$Gene_name[peaks_info$Gene_name == c("RNaseP_nuc")] <- "RNaseP-nuc"

write.table(peaks_info , file="~/Multiomics-ZBF/250106_issaac/gene-peak_pairs_topic8.tsv" , quote = F , sep = "\t")
peaks_info_gc <- StringToGRanges(peaks_info$Peak)

### Choose the peaks  with tbx motif
library(BSgenome.Drerio.UCSC.danRer11)
# library(JASPAR2020)
library(RSQLite)
library(patchwork)
library(TFBSTools)
library(dplyr)
library(Seurat)
library(Signac)
db_path <- "~/JASPAR/JASPAR2024.sqlite3"
pfm <- getMatrixSet(x = db_path , opts = list(collection = "CORE", tax_group = 'vertebrates' , all_versions = FALSE))
#zbf_multi_mesendo <- AddMotifs(zbf_multi_mesendo , genome = BSgenome.Drerio.UCSC.danRer11 ,pfm = pfm)

### Create gene list and peak list splited by chromsome
gr_split <- split(zbf_multi_mesendo@assays$peaks@ranges , seqnames(zbf_multi_mesendo@assays$peaks@ranges))

motif.matrix <- CreateMotifMatrix(
  features = zbf_multi_mesendo@assays$peaks@ranges,
  pwm = pfm,
  genome = BSgenome.Drerio.UCSC.danRer11
)
saveRDS(motif.matrix, file = "~/Multiomics-ZBF/250106_issaac/output/motif.matrix.rds")

motif_names.df <- data.frame(NULL)
for(i in names(pfm))
{
  motif_names.df_tmp <- data.frame(motif_id = pfm[[i]]@ID , motif_name = pfm[[i]]@name)
  motif_names.df <- rbind(motif_names.df , motif_names.df_tmp)
}

tbx_motif <- read.table(file = "~/Multiomics-ZBF/250106_issaac/output/tbx_motif.txt")
motif.matrix_tbx <- motif.matrix[,tbx_motif$V1]
peaks_tbx <- names(which(rowSums(motif.matrix_tbx) == 0))

## SCENT Object
## ATAC matrix
mat_rna <- zbf_multi_mesendo@assays$RNA$counts
mat_atac <- zbf_multi_mesendo@assays$peaks$counts
meta <- zbf_multi_mesendo@meta.data
meta$log_UMI <- log(meta$nCount_RNA)
meta$cell <- meta$cells
##Using the SCENT Object:
SCENT_obj <- CreateSCENTObj(rna = mat_rna, atac = mat_atac, meta.data = meta,
                            peak.info = peaks_info[,c(1,2)],
                            covariates = c("log_UMI"), 
                            celltypes = "stage_cell_type")
## Analysis Topic49 peaks enhancer in Dorsal margin
SCENT_obj_ver1_6hDM <- SCENT_algorithm_opt(object = SCENT_obj, celltype = "6hpf_Dorsal margin", ncores = 6)
write.table(SCENT_obj_ver1_6hDM@SCENT.result , file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_topic49_6hDM.tsv" , sep = "\t" , quote = F)

SCENT_obj_ver1_5hDM <- SCENT_algorithm_opt(object = SCENT_obj, celltype = "5.3hpf_Dorsal margin", ncores = 6)
write.table(SCENT_obj_ver1_5hDM@SCENT.result , file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_topic49_5hDM.tsv" , sep = "\t" , quote = F)

SCENT_obj_ver1_4hDM <- SCENT_algorithm_opt(object = SCENT_obj, celltype = "4.3hpf_Dorsal margin", ncores = 6)
write.table(SCENT_obj_ver1_4hDM@SCENT.result , file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_topic49_4hDM.tsv" , sep = "\t" , quote = F)

# peaks_info[which(!peaks_info$Gene_name %in% rownames(mat_rna)),] %>% View()
# which(gene_ranges$gene_name %in% rownames(zbf_multi_mesendo)) %>% length()
# setdiff(rownames(zbf_multi_mesendo) , gene_ranges$gene_name)

##Example Outputs of the SCENT Object
head(SCENT_obj@rna[1:10,1:2])
head(SCENT_obj@atac[1:10,1:2])
head(SCENT_obj@meta.data)
head(SCENT_obj@peak.info)
str(SCENT_obj)

## SCENT Algorithm: Obtain small list of gene-peak pairs.
#Of the set of peak gene pairs: pick a set of pairs to test: 
#Example: (first 10 gene-peak pairs)
SCENT_obj@peak.info <- SCENT_obj@peak.info[1:10,]
# head(SCENT_obj@peak.info)

## SCENT Algorithm: Options for Regression w/ Bootstrapping.

#Run SCENT algorithm of Tnk cell type and use 6 cores for parallelization:

#Default: Poisson regression and Binarized ATAC counts
# SCENT_obj_ver1 <- SCENT_algorithm(object = SCENT_obj, celltype = "6hpf_Dorsal margin", ncores = 6) 
SCENT_obj_ver1 <- SCENT_algorithm_opt(object = SCENT_obj, celltype = "6hpf_Dorsal margin", ncores = 4) 
# By default settings the above will perform parallelizations using Poisson regression and Binarized counts.

#Option 1: Poisson regression and Non-Binarized ATAC counts
SCENT_obj_ver2 <- SCENT_algorithm(SCENT_obj, "Tnk", 6, regr = "poisson", bin = FALSE)

#Option 2: Negative Binomial regression and Binarized ATAC counts
SCENT_obj_ver3 <- SCENT_algorithm(SCENT_obj, "Tnk", 6, regr = "negbin", bin = TRUE)

#Option 3: Negative Binomial regression and Non-Binarized ATAC counts
SCENT_obj_ver4 <- SCENT_algorithm(SCENT_obj, "Tnk", 6, regr = "negbin", bin = FALSE)


################################### SCENT_algorithm_opt Function
library(data.table)
library(foreach)
library(doParallel)
library(boot)

SCENT_algorithm_opt <- function(object, celltype, ncores, regr = "poisson", bin = TRUE) {
  # 预加载数据到更高效的结构中
  peak_info <- as.data.table(object@peak.info)
  meta_data <- as.data.table(object@meta.data, keep.rownames = "cell")
  setkey(meta_data, cell)
  
  # 预提取ATAC和RNA矩阵（避免重复访问S4对象）
  atac_mat <- as(object@atac, "CsparseMatrix")
  rna_mat <- as(object@rna, "CsparseMatrix")
  cell_names <- colnames(object@atac)
  
  # 预筛选目标细胞类型索引
  target_cells <- meta_data[get(object@celltypes) == celltype, cell]
  if (length(target_cells) == 0) {
    stop("No cells found for celltype: ", celltype)
  }
  # 新增进度提示参数
  total_rows <- nrow(peak_info) 

  # 并行化外层循环
  cl <- makeCluster(ncores, outfile = "")  # 输出重定向到主控制台
  registerDoParallel(cl)
  # registerDoParallel(cores = ncores)
  res <- foreach(n = 1:nrow(peak_info), .combine = rbind, .packages = c("data.table", "boot", "MASS")) %dopar% {
    # 进度提示逻辑
    if (n %% 1 == 0) {
      message(sprintf("Processing %d/%d (%.0f%%)", n, total_rows, n/total_rows*100))
      # sprintf("Processing %d/%d (%.0f%%)", n, total_rows, n/total_rows*100)
    }
    
    gene <- peak_info[n, 1][[1]]
    this_peak <- peak_info[n, 2][[1]]
    
    # 高效提取ATAC数据
    atac_vals <- atac_mat[this_peak, ]
    if (bin) {
      atac_vals <- as.integer(atac_vals > 0)
    }
    atac_dt <- data.table(cell = cell_names, atac = atac_vals)
    
    # 高效提取RNA数据
    exprs_vals <- rna_mat[gene, ]
    expr_dt <- data.table(cell = cell_names, exprs = as.numeric(exprs_vals))
    
    # 快速合并（data.table语法）
    dt <- expr_dt[atac_dt, on = "cell"][meta_data, on = "cell"]
    dt_sub <- dt[cell %in% target_cells]
    
    # 快速计算非零比例
    nonzero_m <- sum(dt_sub$exprs > 0) / nrow(dt_sub)
    nonzero_a <- sum(dt_sub$atac > 0) / nrow(dt_sub)
    
    if (nonzero_m > 0.05 & nonzero_a > 0.05) {
      formula_str <- paste("exprs ~ atac +", paste(object@covariates, collapse = "+"))
      formula <- as.formula(formula_str)
      
      #Estimated Coefficients Obtained without Bootstrapping:
      if(regr == "poisson"){
        base = glm(formula, family = 'poisson', data = dt_sub)
        coefs<-summary(base)$coefficients["atac",]
        assoc <- assoc_poisson
      } else if (regr == "negbin"){
        base = glm.nb(formula, data = dt_sub)
        coefs<-summary(base)$coefficients["atac",]
        assoc <- assoc_negbin
      }
      
      
      # 动态调整bootstrap次数
      bs <- boot(dt_sub, assoc, R = 100, formula = formula, 
                 parallel = "no", ncpus = 1)  # 外层已并行化
      p0 <- basic_p(bs$t0[1], bs$t[,1])
      
      # 调整R_vec和p_vec定义，包含初始R=100
      R_vec <- c(100, 500, 2500, 25000, 50000)
      p_vec <- c(0.1, 0.1, 0.05, 0.01, 0.001)  # R=100对应p<0.1
      
      # 动态调整逻辑（跳过初始的R=100）
      for (i in 2:length(R_vec)) {
        R <- R_vec[i]
        if (p0 < p_vec[i]) {
          bs <- boot(dt_sub, assoc, R = R, formula = formula,
                     parallel = "no", ncpus = 1)
          p0 <- basic_p(bs$t0[1], bs$t[,1])
        } else {
          break
        }
      }
      
      data.frame(gene = gene, peak = this_peak,
                 beta = coefs[1], se = coefs[2], z = coefs[3], p = coefs[4],
                 boot_basic_p = p0)
    } else {
      NULL
    }
  }
  # stopImplicitCluster()
  stopCluster(cl)  # 关闭显式集群（替代 stopImplicitCluster()）
  
  # 保持原始对象结构
  object@SCENT.result <- as.data.frame(res)
  return(object)
}
SCENT_obj_ver1 <- SCENT_algorithm_opt(object = SCENT_obj, celltype = "6hpf_Dorsal margin", ncores = 4) 

CoveragePlot(zbf_multi_mesendo , region = "tbx16l" , group.by = "stage_cell_type" , extend.upstream = 3000 , extend.downstream = 3000)


######## Filtered SCENT result based on FDR < 0.1
input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_raw")
output_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered")

# 获取所有TSV文件列表
tsv_files <- list.files(
  path = input_dir,
  pattern = "\\.tsv$",
  full.names = FALSE
)

# 处理每个文件
for (file_name in tsv_files){
  # 构建完整文件路径
  file_new_name <- sub(pattern = "\\.tsv$" , replacement = "_filtered.tsv" , x = file_name)
  input_path <- file.path(input_dir, file_name)
  output_path <- file.path(output_dir, file_new_name)
  
  # 读取数据
  df <- read.delim(
    file = input_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # 计算FDR校正（假设第3列是p值）
  df$fdr <- p.adjust(df$boot_basic_p, method = "fdr")
  
  # 保存结果
  write.table(
    x = df[df$fdr < 0.1,],
    file = output_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

### Generate uniq peaks
input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered")
output_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_peaks")

tsv_files <- list.files(
  path = input_dir,
  pattern = "\\_filtered.tsv$",
  full.names = FALSE
)

df.new <- data.frame(NULL)
for (file_name in tsv_files){
  # 构建完整文件路径
  file_new_name <- sub(pattern = "\\_filtered.tsv$" , replacement = "_peaks.tsv" , x = file_name)
  input_path <- file.path(input_dir, file_name)
  output_path <- file.path(output_dir, file_new_name)
  
  # 读取数据
  df <- read.delim(
    file = input_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  df$Topic <- (strsplit(x = file_name ,split = "_") %>% unlist)[2]
  df$Topic_celltype <- paste((strsplit(x = file_name ,split = "_") %>% unlist)[2] , (strsplit(x = file_name ,split = "_") %>% unlist)[3] , sep = "_")
  
  df_uniq <- df[!duplicated(df$peak),c(1,2,9,10)]
  df.new <- rbind(df.new , df_uniq)
}
new_df <- df.new %>%
  group_by(peak) %>%
  summarize(
    Topic_celltype = paste(unique(Topic_celltype), collapse = "; "),
    gene = paste(unique(gene), collapse = "; ")
  ) %>%
  ungroup()

peaks_topic47_enhancer.gr <- read.delim(
  file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic47_8hPSM_filtered.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

uniq_topic47_gene <- setdiff(peaks_topic47_enhancer.gr$gene , peaks_topic8_enhancer.gr$gene)
uniq_topic8_gene <- setdiff(peaks_topic8_enhancer.gr$gene , peaks_topic47_enhancer.gr$gene)

markers_PSM <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "8hpf_Presomitic mesoderm (PSM)" , assay = "RNA" , logfc.threshold = 0.5)
markers_5DM <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "5.3hpf_Dorsal margin" , assay = "RNA" , logfc.threshold = 0.5)
markers_5VM <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "5.3hpf_Vertral margin" , assay = "RNA" , logfc.threshold = 0.5)
markers_6DM <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "6hpf_Dorsal margin" , assay = "RNA" , logfc.threshold = 0.5)
markers_6VM <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "6hpf_Vertral margin" , assay = "RNA" , logfc.threshold = 0.5)
markers_END <- FindMarkers(zbf_multi_mesendo, group.by = "stage_cell_type" , ident.1 = "8hpf_Endoderm" , assay = "RNA" , logfc.threshold = 0.5)


markers_topic8_1 <- intersect(rownames(markers_5DM)[markers_5DM$p_val_adj < 0.05] , rownames(markers_5VM)[markers_5VM$p_val_adj < 0.05])
markers_topic8_2 <- intersect(rownames(markers_6DM)[markers_6DM$p_val_adj < 0.05] , rownames(markers_6VM)[markers_6VM$p_val_adj < 0.05])
markers_topic8 <- intersect(markers_topic8_1 , markers_topic8_2)

markers_topic8_withenhancer <- intersect(uniq_topic8_gene , markers_topic8)
markers_topic8_final <- c("fgf3" , "cth1" , "pmepa1" , "mixl1" , "ssbp3a" , "sall3b" , "foxa" , "pitx2" , "dusp4" , "ippk")

# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[1:6] , pt.size = 0 , group.by = "stage_cell_type")
# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[7:12] , pt.size = 0 , group.by = "stage_cell_type")
# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[13:18] , pt.size = 0 , group.by = "stage_cell_type")
# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[19:24] , pt.size = 0 , group.by = "stage_cell_type")
# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[25:30] , pt.size = 0 , group.by = "stage_cell_type")
# VlnPlot(zbf_multi_mesendo , assay = "RNA" , features = markers_topic8_withenhancer[31:35] , pt.size = 0 , group.by = "stage_cell_type")

peaks_topic8_enhancer.gr <- 
  rbind(
  read.delim(
    file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic8_5hDM_filtered.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE),
  read.delim(
    file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic8_6hDM_filtered.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE),
  read.delim(
    file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic8_5hVM_filtered.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE),
  read.delim(
    file = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic8_6hVM_filtered.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE)
  )

# peaks_topic47_enhancer.bed <- read.delim(file = "~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_filtered/SCENT_topic47_8hPSM_filtered.bed",
#                                          sep = "\t",
#                                          header = F,
#                                          stringsAsFactors = FALSE)
# peaks_topic47_enhancer.peaks <- paste(peaks_topic47_enhancer.bed$V1 , peaks_topic47_enhancer.bed$V2 , peaks_topic47_enhancer.bed$V3 , sep = "-")
# peaks_topic47_enhancer.gr <- StringToGRanges(peaks_topic47_enhancer.bed$peak[!duplicated(peaks_topic47_enhancer.bed$peak)])
peaks_topic47_enhancer.gr <- import("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_filtered/SCENT_topic47_8hPSM_filtered.bed")
peaks_topic47_enhancer.gr <- import("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_raw/SCENT_topic47_8hPSM.bed")
export(peaks_topic47_enhancer.gr, con = "~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered/SCENT_topic47_8hPSM_filtered.bed", format = "bed")

CoveragePlot(zbf_multi_mesendo , region = "wnt5b" , group.by = "stage_cell_type" , extend.upstream = 10000, assay = "peaks",extend.downstream = 10000 , 
             region.highlight = StringToGRanges(c("chr1-5537996-5538298" , "chr1-5528995-5529677" , "chr1-5537498-5537942" , "chr1-5542694-5543293")))

peaks_tmp <- df_topic_enhancer[df_topic_enhancer$gene == "sfrp1a",2] %>% StringToGRanges()
peaks_tmp$color <- c(rep("#E64B35" , times = 1) , rep("#00A087" , times = 2))
c2 <- AnnotationPlot(zbf_multi_mesendo , region = "sfrp1a" , extend.upstream = 2000 , extend.downstream = 4000)

c1 <- CoveragePlot(zbf_multi_mesendo , region = c2$plot_env$region , group.by = "stage_cell_type" , 
                   cells = rownames(zbf_multi_mesendo@meta.data)[zbf_multi_mesendo$cell_type %in% c("Dorsal margin" , "Prosterior margin","Vertral margin" , "Presomitic mesoderm (PSM)")],
                   # ranges = peaks_efnb2b)
                   region.highlight = peaks_tmp)

c3 <- BigwigTrack(region = c2$plot_env$region, downsample.rate = 0 , bigwig.scale = "separate",
            bigwig = list(tbx16 = "/mnt/bw/CHIP/zbf_80epi_tbx16_treat_pileup.bw" ,ta = "/mnt/bw/CHIP/zbf_80epi_ta_treat_pileup.bw" , 
                          H3K4me1 = "/mnt/bw/CHIP/zbf_80epi_H3K4me1_treat_pileup.bw",H3K4me3 = "/mnt/bw/CHIP/zbf_80epi_H3K4me3_treat_pileup.bw",
                          H3K27ac = "/mnt/bw/CHIP/zbf_80epi_H3K27ac_treat_pileup.bw") , y_label = "")+NoLegend()
CombineTracks(
  plotlist = list(c3, c1),
  heights = c(1,6)
)

mat <- motif.matrix[c("chr1-5528995-5529677" , "chr1-5537498-5537942" , "chr1-5528507-5528884" , "chr1-5537996-5538298" , "chr1-5542694-5543293"),colSums(motif.matrix[c("chr1-5528995-5529677" , "chr1-5537498-5537942" , "chr1-5528507-5528884" , "chr1-5537996-5538298" , "chr1-5542694-5543293"),]) > 0] %>% as.matrix()
colnames(mat) <- motif_names.df[colnames(mat),2]


CoveragePlot(zbf_multi_mesendo , region = "chrd" , group.by = "stage_cell_type" , extend.upstream = 10000,
             extend.downstream = 10000 , ranges = StringToGRanges(peaks_topic8_enhancer.gr$peak))


CoveragePlot(zbf_multi_mesendo , region = "fn1b" , group.by = "stage_cell_type" , extend.upstream = 10000,
             extend.downstream = 10000 , ranges = StringToGRanges(peaks_topic47_enhancer.gr$peak))


gene <- "fn1b"
c1 <- AnnotationPlot(zbf_multi_mesendo , region = gene , extend.upstream = 20000 , extend.downstream = 20000)

c2 <- BigwigTrack(region = c1$plot_env$region, downsample.rate = 0 , bigwig.scale = "separate",
                  bigwig = list(tbx16 = "/mnt/bw/CHIP/zbf_80epi_tbx16_treat_pileup.bw" ,ta = "/mnt/bw/CHIP/zbf_80epi_ta_treat_pileup.bw" , 
                                H3K4me1 = "/mnt/bw/CHIP/zbf_80epi_H3K4me1_treat_pileup.bw",H3K4me3 = "/mnt/bw/CHIP/zbf_80epi_H3K4me3_treat_pileup.bw",
                                H3K27ac = "/mnt/bw/CHIP/zbf_80epi_H3K27ac_treat_pileup.bw") , y_label = "")+NoLegend()
ranges.show <- subsetByOverlaps(StringToGRanges(peaks_topic8_enhancer.gr$peak), c2$plot_env$region)
c3 <- CoveragePlot(zbf_multi_mesendo , region = c2$plot_env$region , group.by = "stage_cell_type", peaks = F,ranges = ranges.show)

CombineTracks(
  plotlist = list(c2, c3),
  heights = c(1,3)
)

gene <- markers_topic8_final[10]
gene <- "myf5"
c1 <- AnnotationPlot(zbf_multi_mesendo , region = gene , extend.upstream = 20000 , extend.downstream = 20000)
c2 <- BigwigTrack(region = c1$plot_env$region, downsample.rate = 0 , bigwig.scale = "separate",
                  bigwig = list(tbx16 = "/mnt/bw/CHIP/zbf_80epi_tbx16_treat_pileup.bw" ,ta = "/mnt/bw/CHIP/zbf_80epi_ta_treat_pileup.bw" ,
                                H3K4me1 = "/mnt/bw/CHIP/zbf_80epi_H3K4me1_treat_pileup.bw",H3K4me3 = "/mnt/bw/CHIP/zbf_80epi_H3K4me3_treat_pileup.bw",
                                H3K27ac = "/mnt/bw/CHIP/zbf_80epi_H3K27ac_treat_pileup.bw") , y_label = "")+NoLegend()
ranges.show <- subsetByOverlaps(StringToGRanges(peaks_topic8_enhancer.gr$peak), c2$plot_env$region)
c3 <- CoveragePlot(zbf_multi_mesendo , region = c2$plot_env$region , group.by = "stage_cell_type", peaks = F,ranges = ranges.show ,
                   cells = zbf_multi_mesendo$cells[zbf_multi_mesendo$cell_type %in% c("Dorsal margin" , "Prosterior margin" , "Vertral margin" , "Presomitic mesoderm (PSM)")])
# e1 <- ExpressionPlot(zbf_multi_mesendo , features = gene , assay = "RNA")

CombineTracks(
  plotlist = list(c2, c3),
  heights = c(1,3)
)

motif_psm <- FindMarkers(zbf_multi_mesendo_psm , ident.1 = "8hpf_Presomitic mesoderm (PSM)" , assay = "chromvar" , only.pos = T , 
                         mean.fxn = rowMeans, fc.name = "avg_diff" , group.by = "stage_cell_type")

############# Read the peaks in H3K4me1, H3K4me3, H3K27ac
zbf_80epi_H3K4me1 <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_80epi_H3K4me1_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()
zbf_80epi_H3K4me3 <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_80epi_H3K4me3_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()
zbf_80epi_H3K27ac <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_80epi_H3K27ac_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()

zbf_dome_H3K4me1 <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_dome_H3K4me1_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()
zbf_dome_H3K4me3 <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_dome_H3K4me3_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()
zbf_dome_H3K27ac <- with(read.table(file = "~/Multiomics-ZBF/250106_issaac/ChIP_published/zbf_dome_H3K27ac_peaks.broadPeak" , sep = "\t"), 
                          paste(V1 , V2 , V3 , sep = "-")) %>% StringToGRanges()

###### Generate uniq peaks
input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered")
output_dir_1 <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_raw")
output_dir_2 <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_filtered")
tsv_files <- list.files(
  path = input_dir,
  pattern = "\\_filtered.tsv$",
  full.names = FALSE
)

df.new <- data.frame(NULL)
for (file_name in tsv_files){
  # 构建完整文件路径
  file_new_name_1 <- sub(pattern = "\\_filtered.tsv$" , replacement = ".bed" , x = file_name)
  file_new_name_2 <- sub(pattern = "\\_filtered.tsv$" , replacement = "_filtered.bed" , x = file_name)
  input_path <- file.path(input_dir, file_name)
  output_path_1 <- file.path(output_dir_1, file_new_name_1)
  output_path_2 <- file.path(output_dir_2, file_new_name_2)
  # 读取数据
  df <- read.delim(
    file = input_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  # df$Topic <- (strsplit(x = file_name ,split = "_") %>% unlist)[2]
  # df$Topic_celltype <- paste((strsplit(x = file_name ,split = "_") %>% unlist)[2] , (strsplit(x = file_name ,split = "_") %>% unlist)[3] , sep = "_")
  df_uniq <- df$peak[!duplicated(df$peak)] %>% StringToGRanges()
  export(df_uniq, con = output_path_1, format = "bed")
  df_uniq_tmp1 <- df_uniq[intersect(queryHits(findOverlaps(df_uniq , zbf_80epi_H3K4me1)) ,
                                    queryHits(findOverlaps(df_uniq , zbf_80epi_H3K27ac)))]
  df_uniq_tmp11 <- df_uniq_tmp1[!seq_along(df_uniq_tmp1) %in% queryHits(findOverlaps(df_uniq_tmp1 , zbf_80epi_H3K4me3 , minoverlap = 5))]
  
  df_uniq_tmp2 <- df_uniq[intersect(queryHits(findOverlaps(df_uniq , zbf_dome_H3K4me1)) ,
                                    queryHits(findOverlaps(df_uniq , zbf_dome_H3K27ac)))]
  df_uniq_tmp22 <- df_uniq_tmp2[!seq_along(df_uniq_tmp2) %in% queryHits(findOverlaps(df_uniq_tmp2 , zbf_dome_H3K4me3 , minoverlap = 5))]
  
  df_uniq_filtered <- c(df_uniq_tmp11 , df_uniq_tmp22)
  df_uniq_filtered <- df_uniq_filtered[!duplicated(df_uniq_filtered)]
  export(df_uniq_filtered, con = output_path_2, format = "bed")
}

#### Calculate the proportion of candidate enhancer in each topic
### Uniqe peaks in each Topics
# input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_raw")
input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_filtered")
# output_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_Topic_raw")
output_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/Condidate_enhancer_Topic_filtered")
topic_list <- c("topic6","topic7","topic8","topic14","topic47","topic48","topic49")
tsv_files <- list.files(
  path = input_dir,
  pattern = "",
  full.names = FALSE
)
for(i in topic_list)
{
  tsv_files_tmp <- grep(i , tsv_files , value = T)
  df_tmp <- GRanges(NULL)
  for (file_name in tsv_files_tmp){
    bed_tmp <- import(con = paste(input_dir , file_name , sep = "/") , format = "BED" )
    df_tmp <- c(df_tmp , bed_tmp)
  }
  df_tmp <- df_tmp[!duplicated(df_tmp)]
  # export(object = df_tmp , format = "BED" , con = paste(output_dir , "/SCENT_" , i , "_raw.bed" , sep = ""))
  export(object = df_tmp , format = "BED" , con = paste(output_dir , "/SCENT_" , i , "_filtered.bed" , sep = ""))
}

# Import Topic peaks summary
peaks_num <- read.table(file = "~/Multiomics-ZBF/250106_issaac/output/Candidate_enhancer_Topic_summary.txt" , sep = " ")

peaks_num$treat <- factor(rep(x = c("Non-enhancer" , "enhancer") , each = 7) , levels = c("Non-enhancer" , "enhancer"))
peaks_num$topic <- factor(rep(x = c("topic14" , "topic47" , "topic48" , "topic49" , "topic6" , "topic7" , "topic8") , times = 2) , 
                          levels = c("topic14" , "topic48" , "topic49" , "topic8" , "topic6" , "topic7" , "topic47"))
  
peaks_num$value <- peaks_num$V1
peaks_num$value[1:7] <- peaks_num$value[1:7] - peaks_num$value[8:14]

peaks_num$prop <- peaks_num$V1/peaks_num$V1[1:7]
peaks_num$prop[1:7] <- peaks_num$prop[1:7] - peaks_num$prop[8:14]

ggplot(peaks_num, aes(x = topic, y = prop , fill = treat)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Topic", x = "Topic", y = "Peaks_Counts")+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))+
  theme_minimal()+
  scale_fill_manual(values = c("Non-enhancer" = "#4E79A7", "enhancer" = "#E15759"))+
  geom_text(
    aes(label = sprintf("%.1f%%", prop*100)),  # 将比例转换为百分比格式
    position = position_stack(vjust = 0.5),    # 标签居中显示在堆叠区域
    color = "white",                           # 设置高对比度文字颜色
    size = 3.5,                                # 适当调整字号
    fontface = "bold"                          # 加粗字体提升可读性
  )
  
ggplot(peaks_num, aes(x = topic, y = value , fill = treat)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Topic", x = "Topic", y = "Peaks_Counts")+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))+
  theme_minimal()+
  scale_fill_manual(values = c("Non-enhancer" = "#4E79A7", "enhancer" = "#E15759"))+
  geom_text(
    aes(label = value),  # 将比例转换为百分比格式
    position = position_stack(vjust = 0.5),    # 标签居中显示在堆叠区域
    color = "white",                           # 设置高对比度文字颜色
    size = 3.5,                                # 适当调整字号
    fontface = "bold"                          # 加粗字体提升可读性
  )

##### GO analysis of candidate enhancers' target genes
## Generate target_gene list and df
library(clusterProfiler)
library(org.Dr.eg.db)

input_dir <- path.expand("~/Multiomics-ZBF/250106_issaac/output/SCENT_filtered")
tsv_files <- list.files(
  path = input_dir,
  pattern = "",
  full.names = FALSE
)
target_gene_list <- list()
target_gene_df <- data.frame(NULL)
for(file_name in tsv_files)
{
  list_name <- sub("^SCENT_(.*)_filtered\\..*$", "\\1", file_name)
  df <- read.delim(
    file = paste(input_dir , file_name , sep = "/"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  df <- df[df$beta > 0,]
  target_gene_tmp <- df$gene[!duplicated(df$gene)]
  if(length(target_gene_tmp) > 0)
  {
    target_gene_list[[list_name]] <- target_gene_tmp
    target_gene_df <- rbind(target_gene_df , data.frame(Gene_name = target_gene_tmp, topic_celltype = list_name))
  }else{
    next
  }
}

### Venn plot of each topics
library(ggVennDiagram)
ggVennDiagram(
  list(A = target_gene_list[['topic8_4hDM']], B = target_gene_list[['topic8_5hDM']], C = target_gene_list[['topic8_6hDM']] , 
       D = target_gene_list[['topic8_5hVM']], E = target_gene_list[['topic8_6hVM']]),  # 输入数据为命名列表
  label = "both",             # 显示数量和百分比
  label_alpha = 0.5,          # 标签背景透明度
  category.names = c("topic8_4hDM", "topic8_5hDM", "topic8_6hDM" , "topic8_5hVM", "topic8_6hVM")  # 自定义分组名称
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#2196F3") +  # 颜色渐变
  theme(legend.position = "none")  # 隐藏图例

ggVennDiagram(
  list(A = target_gene_list[['topic8_4hDM']], B = target_gene_list[['topic8_5hDM']], C = target_gene_list[['topic8_6hDM']]),  # 输入数据为命名列表
  label = "both",             # 显示数量和百分比
  label_alpha = 0.5,          # 标签背景透明度
  category.names = c("topic8_4hDM", "topic8_5hDM", "topic8_6hDM")  # 自定义分组名称
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#2196F3") +  # 颜色渐变
  theme(legend.position = "none")  # 隐藏图例

ggVennDiagram(
  list(B = target_gene_list[['topic8_5hVM']], C = target_gene_list[['topic8_6hVM']]),  # 输入数据为命名列表
  label = "both",             # 显示数量和百分比
  label_alpha = 0.5,          # 标签背景透明度
  category.names = c("topic8_5hVM", "topic8_6hVM")  # 自定义分组名称
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#2196F3") +  # 颜色渐变
  theme(legend.position = "none")  # 隐藏图例

## Perform GO
background_genes <- zbf_multi_mesendo@assays$RNA$counts[which(rowSums(zbf_multi_mesendo@assays$RNA$counts) > 0),] %>% rownames()

table(target_gene_df$topic_celltype) %>% names()
topic_celltype <- "topic47_8hPSM"
target_gene <- target_gene_df$Gene_name[target_gene_df$topic_celltype == topic_celltype]
# target_gene <- target_gene_df$Gene_name[grep("topic8" , target_gene_df$topic_celltype)]
# target_gene <- target_gene[!duplicated(target_gene)]

go_enrich <- enrichGO(
  gene = target_gene, 
  universe = background_genes, 
  OrgDb = org.Dr.eg.db,       
  keyType = "SYMBOL",         
  ont = "ALL",                # 本体论：BP（生物学过程）、MF（分子功能）、CC（细胞组分）、ALL
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
# dotplot(go_enrich, showCategory = 20, title = paste("GO" , topic_celltype , sep = "_"))
barplot(go_enrich, showCategory = 20, title = paste("GO" , topic_celltype , sep = "_"))

target_gene_entrez <- bitr(target_gene , fromType = "SYMBOL" , toType = "ENTREZID" , OrgDb = org.Dr.eg.db)
kegg_enrich <- enrichKEGG(
  gene = target_gene_entrez$ENTREZID, 
  organism = "dre" , 
  pvalueCutoff = 0.05 , 
  pAdjustMethod = "BH")
# dotplot(kegg_enrich, showCategory = 20, title = paste("KEGG" , topic_celltype , sep = "_"))
# barplot(kegg_enrich, showCategory = 20, title = paste("KEGG" , topic_celltype , sep = "_"))

# topic8_go_list <- list()
# topic8_kegg_list <- list()
topic8_go_list[[topic_celltype]] <- go_enrich@result$Description[go_enrich@result$p.adjust < 0.05]
topic8_kegg_list[[topic_celltype]] <- kegg_enrich@result$Description[kegg_enrich@result$p.adjust < 0.05]

ggVennDiagram(
  list(B = topic8_go_list[['topic8_5hDM']], C = topic8_go_list[['topic8_6hDM']] , 
       D = topic8_go_list[['topic8_5hVM']], E = topic8_go_list[['topic8_6hVM']]),  # 输入数据为命名列表
  label = "both",             # 显示数量和百分比
  label_alpha = 0.5,          # 标签背景透明度
  category.names = c("topic8_5hDM", "topic8_6hDM" , "topic8_5hVM", "topic8_6hVM")  # 自定义分组名称
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#2196F3") +  # 颜色渐变
  theme(legend.position = "none")  # 隐藏图例

pathway_1 <- intersect(x = topic8_go_list$topic8_5hDM , y = topic8_go_list$topic8_6hDM)
pathway_2 <- intersect(x = topic8_go_list$topic8_5hVM , y = topic8_go_list$topic8_6hVM)
pathway <- intersect(x = pathway_1 , y = pathway_2)

combined <- unlist(topic8_go_list)
freq <- table(combined)
freq <- freq[freq >= 4]

# Combine common pathway
go_enrich_topic8_common <- data.frame(NULL)

for(topic_celltype in c("topic8_5hDM" , "topic8_5hVM" , "topic8_6hDM" , "topic8_6hVM"))
{
  celltype <- substr(topic_celltype , start = 8 , stop = 11)
  target_gene <- target_gene_df$Gene_name[target_gene_df$topic_celltype == topic_celltype]
  # target_gene <- target_gene_df$Gene_name[grep("topic8" , target_gene_df$topic_celltype)]
  # target_gene <- target_gene[!duplicated(target_gene)]
  
  go_enrich <- enrichGO(
    gene = target_gene, 
    universe = background_genes, 
    OrgDb = org.Dr.eg.db,       
    keyType = "SYMBOL",         
    ont = "ALL",                # 本体论：BP（生物学过程）、MF（分子功能）、CC（细胞组分）、ALL
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  go_enrich_topic8_tmp <- go_enrich@result[go_enrich@result$Description %in% names(freq),]
  go_enrich_topic8_tmp$celltype <- celltype
  go_enrich_topic8_common <- rbind(go_enrich_topic8_common , go_enrich_topic8_tmp)
}

go_enrich_topic8_common$pathway <- paste(go_enrich_topic8_common$celltype , go_enrich_topic8_common$Description , sep = "_")

### Summary genes both in topic47 and topic 8
target_gene_df_topic8 <- target_gene_df[target_gene_df$topic_celltype %in% c("topic8_4hDM","topic8_5hDM","topic8_5hVM","topic8_6hDM","topic8_6hVM"),]
target_gene_df_topic47 <- target_gene_df[target_gene_df$topic_celltype == "topic47_8hPSM",]

table(target_gene_df_topic8$Gene_name)[table(target_gene_df_topic8$Gene_name) > 4] %>% names() %in% target_gene_df_topic47$Gene_name



