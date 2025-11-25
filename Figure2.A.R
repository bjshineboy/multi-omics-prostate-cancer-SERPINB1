# ====================== Optimized Survival Analysis ======================
# Objective: Analyze the relationship between SERPINB1 methylation level and 
#            progression-free survival (PFS) in prostate cancer
# =========================================================================

# Load required R packages
# =========================================================================
suppressPackageStartupMessages({
  library(data.table)    # Efficient reading of large datasets
  library(dplyr)         # Data manipulation
  library(tidyr)         # Data tidying
  library(survival)      # Survival analysis
  library(survminer)     # Survival curve visualization
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Methylation array annotation
})


setwd("Figure2_Data")

# 1. Load methylation data
# =========================================================================
cat("\n===== Step 1: Load methylation data =====\n")
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

# Efficiently read compressed data using fread
methylation_data <- tryCatch({
  data.table::fread(
    cmd = paste("zcat", meth_file),
    header = TRUE,
    sep = "\t",
    data.table = FALSE,
    showProgress = FALSE  # Hide progress bar
  )
}, error = function(e) {
  stop("Failed to load methylation data: ", e$message)
})

# Validate basic data structure
if (nrow(methylation_data) < 2 || ncol(methylation_data) < 2) {
  stop("Invalid methylation data format: rows=", nrow(methylation_data), " columns=", ncol(methylation_data))
}

# Extract sample IDs and CpG site names
sample_names <- colnames(methylation_data)[-1] 
cpg_names <- as.character(methylation_data[-1, 1])  

# Create Beta value matrix
beta_matrix <- as.matrix(methylation_data[-1, -1])
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names

# Convert to data frame with CpG column
beta_df <- as.data.frame(beta_matrix) %>%
  tibble::rownames_to_column(var = "CpG")

cat("✅ Methylation data loaded successfully\n")
cat("   - Samples:", ncol(beta_matrix), "\n")
cat("   - CpG sites:", nrow(beta_matrix), "\n")
cat("   - Beta value range:", round(range(as.vector(beta_matrix), na.rm = TRUE), 3), "\n")

# 2. Load clinical data
# =========================================================================
cat("\n===== Step 2: Load clinical data =====\n")
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"
surv_file <- "data/survival_PRAD_survival.txt"

# Key clinical variables
# Read clinical data (only essential columns to conserve memory)
clinical_data <- fread(
  clinical_file,
  select = c("sampleID",   "sample_type"),
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")  # Consistent NA markers
)

surv_file_data <- fread(
  surv_file,
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")  # Consistent NA markers
)


surv_file_data_merge_clin_tumor <- merge(surv_file_data, clinical_data, by.x="sample", by.y ="sampleID", all=TRUE )%>%
  filter(sample_type== "Primary Tumor")  


# 3. Identify SERPINB1-related CpG sites
# =========================================================================
cat("\n===== Step 3: Identify SERPINB1-related CpG sites =====\n")
# Load 450k array annotation data
anno_data <- IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19

# Extract annotation information
anno_df <- as.data.frame(getAnnotation(anno_data))

# Filter CpG sites related to SERPINB1 gene
cpgs_serpinb1 <- anno_df %>%
  filter(grepl("SERPINB1\\b", UCSC_RefGene_Name)) %>%
  dplyr::select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_Island) %>%
  filter(UCSC_RefGene_Group != "Body")

cpg_list <- unique(cpgs_serpinb1$Name)

cat("✅ Identified", length(cpg_list), "SERPINB1-related CpG sites\n")
cat("   - Related genomic regions:", paste(unique(cpgs_serpinb1$UCSC_RefGene_Group), collapse = ", "), "\n")

# 4. Calculate average SERPINB1 methylation level per sample
# =========================================================================
cat("\n===== Step 4: Calculate average methylation level =====\n")

# Extract methylation data for SERPINB1-related CpGs
meth_serpinb1 <- beta_df %>%
  filter(CpG %in% cpg_list)

if (nrow(meth_serpinb1) == 0) {
  stop("No matching SERPINB1-related CpG sites found")
}

# Convert to long format
meth_long <- meth_serpinb1 %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  )

# Calculate average methylation per sample
avg_meth <- meth_long %>%
  group_by(sample_id) %>%
  summarise(
    Avg_Beta = mean(Beta_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Extract patient ID (first 12 characters of sample ID)
  mutate(patient_id = substr(sample_id, 1, 12))

cat("✅ Average methylation calculation complete\n")
cat("   - Samples with methylation data:", nrow(avg_meth), "\n")

# 5. Data integration and grouping
# =========================================================================
cat("\n===== Step 5: Integrate data and create groups =====\n")

surv_file_data_merge_clin_tumor_beta_group <- surv_file_data_merge_clin_tumor %>%
  left_join(avg_meth %>% dplyr::select(sample_id, Avg_Beta), by = c("sample" = "sample_id")) %>%
  # Create methylation groups based on median
  mutate(
    Methylation_Group = case_when(
      Avg_Beta > 0.75 ~ "High Methylation",
      Avg_Beta <= 0.75 ~ "Low Methylation",
      TRUE ~ NA_character_  # Will be filtered out
    ),
    Methylation_Group = factor(Methylation_Group, levels = c("Low Methylation", "High Methylation"))
  )
 
# Group counts
group_counts <- surv_file_data_merge_clin_tumor_beta_group %>%
  dplyr::count(Methylation_Group)

cat("✅ Data integration complete\n")
cat("   - Final analysis cohort:", nrow(surv_file_data_merge_clin_tumor_beta_group), "\n")
cat("   - Group distribution:\n")
for (i in 1:nrow(group_counts)) {
  cat("     -", group_counts$Methylation_Group[i], ":", group_counts$n[i], "samples\n")
}

# 6. Survival analysis
# =========================================================================
cat("\n===== Step 6: Survival analysis =====\n")
library(survival)
library(survminer)

# 加载必要包
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)
 

# 创建生存对象列表
fit_list <- list(
  OS = survfit(Surv(OS.time/365, OS) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group),
  DSS = survfit(Surv(DSS.time/365, DSS) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group),
  PFI = survfit(Surv(PFI.time/365, PFI) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group)
)


# 专业医学主题
medical_theme <- theme(
  text = element_text(family = "Arial"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 10),
  legend.position = "bottom",
  #panel.grid.major = element_line(color = "gray90", size = 0.2),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(1, 1, 1, 1), "cm")
)

# 专业配色方案
clinical_palette <- c("#1B9E77", "#D95F02")  # 绿色-橙色组合

# 创建生存曲线绘制函数
create_survplot <- function(fit, title) {
  ggsurvplot(fit,
             data = surv_file_data_merge_clin_tumor_beta_group,
             palette = clinical_palette,
             title = title,
             xlab = "Time (Years)",
             ylab = "Survival Probability",
             pval = TRUE,
             pval.size = 4,
             pval.coord = c(1, 0.15),
             risk.table = TRUE,
             risk.table.height = 0.25,
             risk.table.title = "",
             risk.table.y.text = FALSE,
             risk.table.fontsize = 3.5,
             tables.theme = theme_cleantable(),
             ncensor.plot = FALSE,
             ggtheme = medical_theme,
             break.time.by = 1,
             legend.labs = levels(surv_file_data_merge_clin_tumor_beta_group$Methylation_Group),
             legend.title = "",
             surv.plot.height = 0.7,
             tables.height = 0.3)
}



# 生成各终点图形
plots <- list(
  OS = create_survplot(fit_list$OS, "Overall Survival"),
  DSS = create_survplot(fit_list$DSS, "Disease-Specific Survival"),
  PFI = create_survplot(fit_list$PFI, "Progression-Free Interval")
)

# 整合图形布局
final_plot <- wrap_plots(
  plots$OS$plot + theme(legend.position = "none"),
  plots$DSS$plot + theme(legend.position = "none"),
  plots$PFI$plot,
  plots$OS$table,
  plots$DSS$table,
  plots$PFI$table,
  design = "
  ABC
  DEF
  ",
  widths = c(1, 1, 1),  # 等宽分布
  heights = c(1.5, 0.4)  # 主图与风险表高度比
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(10, 0, 10, 0),
    plot.margin = unit(c(0.5, 0.3, 0.5, 0.3), "cm"),  # 减少左右边距
    legend.text = element_text(size = 10)
  )

# 添加全局标题和注释
final_plot <- final_plot +
  plot_annotation(
    title = "Triple Endpoint Survival Analysis",
    subtitle = "Parallel Comparison of Methylation Effects",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# 确保y轴标签对齐
final_plot <- final_plot &
  theme(
    axis.title.y = element_text(vjust = 0.5, margin = margin(r = 10)),
    plot.title.position = "plot"
  )

# 生成统计汇总表
results_summary <- map_dfr(names(fit_list), function(x) {
  s <- summary(fit_list[[x]])
  data.frame(
    Endpoint = x,
    Group = levels(surv_file_data_merge_clin_tumor_beta_group$Methylation_Group),
    Median_Survival = s$table[, "median"],
    Year1_Survival = s$surv[which.min(abs(s$time - 1))],
    Year3_Survival = s$surv[which.min(abs(s$time - 3))],
    P_value = surv_pvalue(fit_list[[x]])$pval
  )
})

# 打印统计结果
print(results_summary)

# PDF输出（矢量图）
pdf("TCGA_PRAD_Methylation_Survival2.pdf", width = 7500, height = 3400)
print(final_plot)
dev.off()

# PNG输出（高分辨率）
png("TCGA_PRAD_Methylation_Survival.png", width = 7500, height = 3400, res = 300)
print(final_plot)
dev.off()

write.csv(surv_file_data_merge_clin_tumor_beta_group, "surv_file_data_merge_clin_tumor_beta_group.csv",  row.names = FALSE)


###########PFI
# 生存分析
pfi_survdiff <- survdiff(Surv(PFI.time, PFI) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group)

# Cox回归
pfi_cox <- coxph(Surv(PFI.time, PFI) ~ Methylation_Group , 
      data = surv_file_data_merge_clin_tumor_beta_group)


pfi_survdiff
pfi_cox


# 中位生存时间
print(fit_list$PFI)

# HR与P值
summary(pfi_cox)


###############OS

# 生存分析
pfi_survdiff_OS <- survdiff(Surv(OS.time, OS) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group)

# Cox回归
pfi_cox_OS <- coxph(Surv(OS.time, OS) ~ Methylation_Group , 
                 data = surv_file_data_merge_clin_tumor_beta_group)


pfi_survdiff_OS
pfi_cox_OS


# 中位生存时间
print(fit_list$OS)

# HR与P值
summary(pfi_cox_OS)

###########################DSS
###############DSS

# 生存分析
pfi_survdiff_DSS <- survdiff(Surv(DSS.time, DSS) ~ Methylation_Group, data = surv_file_data_merge_clin_tumor_beta_group)

# Cox回归
pfi_cox_DSS <- coxph(Surv(DSS.time, DSS) ~ Methylation_Group , 
                    data = surv_file_data_merge_clin_tumor_beta_group)


pfi_survdiff_DSS
pfi_cox_DSS


# 中位生存时间
print(fit_list$DSS)

# HR与P值
summary(pfi_cox_DSS)


