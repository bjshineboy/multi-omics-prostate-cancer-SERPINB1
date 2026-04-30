# ==============================================================================
# Prostate Cancer Methylation Data Analysis: SERPINB1 Multi-Analysis Figure
# 
# Objective: Generate a combined figure with three subplots:
#   A: Correlation between SERPINB1-related CpG sites and Gleason Score (bar plot)
#   B: Correlation between SERPINB1 methylation and E-MAPK signaling genes (bar plot)
#   C: Correlation between SERPINB1 promoter methylation and expression (scatter plot)
#
# Author: Bioinformatics Analysis Team
# Date: 2025-06-27
# ==============================================================================

# --------------------- 1. Load Required Packages ---------------------
library(data.table)     # Efficient reading of large data files
library(dplyr)          # Data manipulation and piping
library(tidyr)          # Data transformation
library(ggplot2)        # Basic visualization
library(cowplot)        # Plot combining

setwd("/work/singleCell/prostate_singecell_20250627/GSE_Raw/sci/SERPNB1_sci_20250829/code/2/Figure3_Data")

# ===================== PART 1: Gleason Score Correlation =====================
cat("PART 1: SERPINB1 CpG Methylation vs Gleason Score Correlation\n")

# --------------------- Read and Preprocess Methylation Data ---------------------
meth_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

cat("Step 1.1: Reading methylation data...\n")
methylation_data <- tryCatch(
  {
    fread(
      cmd = paste("zcat", meth_file),
      header = TRUE,
      sep = "\t",
      data.table = FALSE,
      showProgress = TRUE
    )
  },
  error = function(e) {
    stop(paste("Failed to read methylation data:", e$message))
  }
)

# Process methylation matrix
cpg_names <- as.character(methylation_data[-1, 1])
sample_names <- colnames(methylation_data)[-1]
beta_matrix <- as.matrix(methylation_data[-1, -1])
mode(beta_matrix) <- "numeric"
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names
beta_df <- as.data.frame(beta_matrix) %>% rownames_to_column(var = "CpG")

cat("Methylation data loaded: ", nrow(beta_matrix), "CpG sites, ", ncol(beta_matrix), "samples\n")

# --------------------- Read Clinical Data ---------------------
clinical_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

cat("Step 1.2: Reading clinical data...\n")
clinical_data <- tryCatch(
  {
    fread(
      clinical_file,
      header = TRUE,
      sep = "\t",
      data.table = FALSE,
      na.strings = c("", "NA", "[Not Available]")
    )
  },
  error = function(e) {
    stop(paste("Failed to read clinical data:", e$message))
  }
)

clinical_df <- clinical_data %>%
  select(
    sample_id = "sampleID",
    Gleason = "gleason_score",
    sample_type
  ) %>%
  mutate(
    Gleason = as.numeric(as.character(Gleason)),
  ) %>%
  filter(
    !is.na(Gleason),
    sample_type == "Primary Tumor",
    Gleason >= 5, Gleason <= 10
  )

cat("Clinical data loaded: ", nrow(clinical_df), "primary tumor samples\n")

# --------------------- Identify SERPINB1 CpG Sites ---------------------
cat("Step 1.3: Identifying SERPINB1-related CpG sites...\n")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # 450K array annotation
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cpgs_serpinb1 <- as.data.frame(anno) %>%
  filter(grepl("SERPINB1\\b", UCSC_RefGene_Name))

cpg_list <- cpgs_serpinb1$Name
cat("Found", length(cpg_list), "CpG sites related to SERPINB1\n")

# Extract CpG site details
cpgs_serpinb1_info <- cpgs_serpinb1 %>%
  select(
    Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group
  ) %>%
  mutate(
    Position = paste0(chr, ":", pos),
    Gene_Region = case_when(
      grepl("TSS", UCSC_RefGene_Group) ~ "Promoter",
      grepl("Body", UCSC_RefGene_Group) ~ "Gene Body",
      grepl("5'UTR", UCSC_RefGene_Group) ~ "5'UTR",
      grepl("3'UTR", UCSC_RefGene_Group) ~ "3'UTR",
      TRUE ~ "Other"
    )
  )

# --------------------- Integrate Data ---------------------
cat("Step 1.4: Integrating methylation and clinical data...\n")
meth_serpinb1 <- beta_df %>%
  filter(CpG %in% cpg_list) %>%
  select(CpG, all_of(clinical_df$sample_id)) %>%
  filter(rowSums(!is.na(select(., -CpG))) > 0)

meth_long <- meth_serpinb1 %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  filter(!is.na(Beta_value))

# Merge with clinical data
meth_long <- meth_long %>%
  mutate(patient_id = substr(sample_id, 1, 12))

clinical_df <- clinical_df %>%
  mutate(patient_id = substr(sample_id, 1, 12))

merged_data_p1 <- meth_long %>%
  inner_join(
    clinical_df %>% select(patient_id, Gleason, sample_type),
    by = "patient_id"
  ) %>%
  left_join(
    cpgs_serpinb1_info,
    by = c("CpG" = "Name")
  )

cat("Integrated data: ", nrow(merged_data_p1), "rows, ", 
    n_distinct(merged_data_p1$patient_id), "patients, ",
    n_distinct(merged_data_p1$CpG), "CpG sites\n")

# --------------------- Calculate Correlations ---------------------
cat("Step 1.5: Calculating Spearman correlations...\n")
cor_results_p1 <- merged_data_p1 %>%
  group_by(CpG, Position, Gene_Region) %>%
  summarise(
    Spearman_corr = cor.test(Beta_value, Gleason, method = "spearman")$estimate,
    P_value = cor.test(Beta_value, Gleason, method = "spearman")$p.value,
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Direction = ifelse(Spearman_corr > 0, "Positive", "Negative"),
    CpG_label = paste0(CpG, " (", Gene_Region, ")")
  ) %>%
  arrange(desc(abs(Spearman_corr)))

write.csv(cor_results_p1, "SERPINB1_Gleason_correlation_results.csv", row.names = FALSE)
cat("Saved correlation results to: SERPINB1_Gleason_correlation_results.csv\n")

# ===================== PART 2: E-MAPK Gene Correlation =====================
cat("PART 2: SERPINB1 Methylation vs E-MAPK Gene Expression Correlation\n")

# --------------------- Read Expression Data ---------------------
expr_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_HiSeqV2.gz"

cat("Step 2.1: Reading expression data...\n")
expr_data <- fread(
  cmd = paste("zcat", expr_file),
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

# Process expression matrix
gene_names_expr <- as.character(expr_data[, 1])
sample_names_expr <- colnames(expr_data)[-1]
expr_matrix <- as.matrix(expr_data[, -1])
rownames(expr_matrix) <- gene_names_expr
colnames(expr_matrix) <- sample_names_expr
expr_df <- as.data.frame(expr_matrix) %>% rownames_to_column(var = "Gene")

cat("Expression data loaded: ", nrow(expr_matrix), "genes, ", ncol(expr_matrix), "samples\n")

# --------------------- Calculate Average SERPINB1 Methylation ---------------------
cat("Step 2.2: Calculating average SERPINB1 methylation...\n")
avg_meth <- beta_df %>%
  filter(CpG %in% cpg_list) %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    SERPINB1_Beta = mean(Beta_value, na.rm = TRUE)
  )

# --------------------- Extract E-MAPK Genes ---------------------
cat("Step 2.3: Extracting E-MAPK signaling axis genes...\n")
target_genes <- c("MMP2", "MMP9", "SNAI1", "TWIST1", "MAPK3", "MAPK1")

expr_target <- expr_df %>%
  filter(Gene %in% target_genes) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample_id",
    values_to = "Expression"
  ) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Expression
  )

# --------------------- Integrate and Calculate Correlations ---------------------
cat("Step 2.4: Calculating correlations with E-MAPK genes...\n")
common_samples <- intersect(avg_meth$sample_id, expr_target$sample_id)
avg_meth_common <- avg_meth %>% filter(sample_id %in% common_samples)
expr_common <- expr_target %>% filter(sample_id %in% common_samples)

merged_data_p2 <- inner_join(avg_meth_common, expr_common, by = "sample_id")

# Calculate correlations
cor_results_p2 <- data.frame()
for (gene in target_genes) {
  if (gene %in% colnames(merged_data_p2)) {
    cor_test <- cor.test(merged_data_p2$SERPINB1_Beta, merged_data_p2[[gene]], 
                         method = "spearman", use = "complete.obs")
    cor_results_p2 <- rbind(cor_results_p2, data.frame(
      Gene = gene,
      Spearman_corr = cor_test$estimate,
      P_value = cor_test$p.value,
      Sample_size = sum(complete.cases(merged_data_p2$SERPINB1_Beta, merged_data_p2[[gene]]))
    ))
  }
}

cor_results_p2 <- cor_results_p2 %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(desc(abs(Spearman_corr)))

write.csv(cor_results_p2, "SERPINB1_EMAPK_correlation_results.csv", row.names = FALSE)
cat("Saved E-MAPK correlation results to: SERPINB1_EMAPK_correlation_results.csv\n")

# ===================== PART 3: Methylation-Expression Correlation =====================
cat("PART 3: SERPINB1 Promoter Methylation vs Expression Correlation\n")

# --------------------- Extract Promoter CpG Sites ---------------------
cat("Step 3.1: Extracting SERPINB1 promoter CpG sites...\n")
cpgs_promoter_serpinb1 <- as.data.frame(anno) %>%
  filter(
    grepl("SERPINB1\\b", UCSC_RefGene_Name),
    grepl("TSS1500|TSS200|1stExon|5'UTR", UCSC_RefGene_Group)
  )

cpg_list_promoter <- unique(cpgs_promoter_serpinb1$Name)
cat("Found", length(cpg_list_promoter), "promoter CpG sites for SERPINB1\n")

# Calculate average promoter methylation
avg_meth_promoter <- beta_df %>%
  filter(CpG %in% cpg_list_promoter) %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    SERPINB1_Promoter_Beta = mean(Beta_value, na.rm = TRUE),
    CpG_count = n()
  )

# Extract SERPINB1 expression
expr_serpinb1 <- expr_df %>%
  filter(Gene == "SERPINB1") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample_id",
    values_to = "SERPINB1_Expression"
  ) %>%
  select(-Gene)

# --------------------- Classify Samples ---------------------
cat("Step 3.2: Classifying samples as Cancer/Normal...\n")
clinical_classify <- clinical_data %>%
  select(sample_id = "sampleID", sample_type) %>%
  mutate(
    Group = case_when(
      grepl("tumor|primary", sample_type, ignore.case = TRUE) ~ "Cancer",
      grepl("normal|solid tissue", sample_type, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Group %in% c("Cancer", "Normal"))

# --------------------- Integrate Data ---------------------
cat("Step 3.3: Integrating methylation, expression, and clinical data...\n")
merged_data_p3 <- avg_meth_promoter %>%
  inner_join(expr_serpinb1, by = "sample_id") %>%
  inner_join(clinical_classify, by = "sample_id")

# Calculate correlations by group
cor_results_p3 <- merged_data_p3 %>%
  group_by(Group) %>%
  summarise(
    Spearman_corr = cor(SERPINB1_Promoter_Beta, SERPINB1_Expression, 
                        method = "spearman", use = "complete.obs"),
    P_value = cor.test(SERPINB1_Promoter_Beta, SERPINB1_Expression, 
                       method = "spearman", exact = FALSE)$p.value,
    N = n(),
    .groups = "drop"
  )

write.csv(cor_results_p3, "SERPINB1_Meth_Expr_correlation_results.csv", row.names = FALSE)
cat("Saved methylation-expression correlation results\n")

# ===================== CREATE COMBINED FIGURE =====================
cat("Creating combined figure with three subplots (A, B, C)\n")

# --------------------- Create Plot A: Bar Plot (CpG vs Gleason) ---------------------
cat("Creating Plot A: Gleason correlation bar plot...\n")


n_cpgs_to_show <- nrow(cor_results_p1)
if (n_cpgs_to_show > 20) {
  #  
  n_cpgs_to_show <- 20
  cor_results_p1_display <- cor_results_p1 %>%
    arrange(desc(abs(Spearman_corr))) %>%
    head(n_cpgs_to_show)
} else {
  cor_results_p1_display <- cor_results_p1
}

#  
gene_region_colors <- c(
  "Promoter" = "#FF6B6B", 
  "Gene Body" = "#4ECDC4", 
  "5'UTR" = "#45B7D1", 
  "3'UTR" = "#96CEB4",
  "Other" = "#FFEAA7"
)

#  
p_a <- ggplot(cor_results_p1_display, aes(x = reorder(CpG_label, Spearman_corr), 
                                          y = Spearman_corr, 
                                          fill = Gene_Region)) +
  geom_bar(stat = "identity", width = 0.7) +
  # 
  geom_text(aes(label = sprintf("%.3f", Spearman_corr), 
                y = ifelse(Spearman_corr > 0, Spearman_corr + 0.02, Spearman_corr - 0.02)),
            size = 2.5, hjust = ifelse(cor_results_p1_display$Spearman_corr > 0, -0.2, 1.2)) +
  scale_fill_manual(values = gene_region_colors) +
  labs(
    title = "SERPINB1 CpG Methylation\nvs Gleason Score",
    x = "CpG Site (Genomic Region)",
    y = "Spearman Correlation (ρ)",
    fill = "Gene Region"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "#2C3E50"),
    axis.title = element_text(size = 9, face = "bold", color = "#2C3E50"),
    axis.text = element_text(size = 7, color = "#34495E"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_line(color = "#ECF0F1"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5, "mm")
  ) +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3)

# --------------------- Create Plot B: Bar Plot (E-MAPK Genes) ---------------------
cat("Creating Plot B: E-MAPK correlation bar plot...\n")

# 创建B图：柱状图
p_b <- ggplot(cor_results_p2, aes(x = reorder(Gene, Spearman_corr), 
                                  y = Spearman_corr, 
                                  fill = ifelse(Spearman_corr > 0, "Positive", "Negative"))) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Significance, 
                y = ifelse(Spearman_corr > 0, Spearman_corr + 0.05, Spearman_corr - 0.05)),
            size = 3, vjust = ifelse(cor_results_p2$Spearman_corr > 0, -0.5, 1.2)) +
  scale_fill_manual(values = c("Positive" = "#A23B72", "Negative" = "#2E86AB")) +
  labs(
    title = "E-MAPK Signal Axis Genes\nvs SERPINB1 Methylation",
    x = "Gene",
    y = "Spearman Correlation (ρ)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "#2C3E50"),
    axis.title = element_text(size = 9, face = "bold", color = "#2C3E50"),
    axis.text = element_text(size = 8, color = "#34495E"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major = element_line(color = "#ECF0F1"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5, "mm")
  ) +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3)

# --------------------- Create Plot C: Scatter Plot ---------------------
cat("Creating Plot C: Methylation-expression scatter plot...\n")

# 创建散点图
p_c <- ggplot(merged_data_p3, 
              aes(x = SERPINB1_Promoter_Beta, 
                  y = SERPINB1_Expression, 
                  color = Group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.5) +
  facet_wrap(~ Group, scales = "free") +
  labs(
    title = "SERPINB1 Promoter Methylation\nvs Gene Expression",
    x = "Promoter Methylation (β-value)",
    y = "SERPINB1 Expression"
  ) +
  scale_color_manual(
    values = c("Cancer" = "#E15759", "Normal" = "#4E79A7")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "#2C3E50"),
    axis.title = element_text(size = 9, face = "bold", color = "#2C3E50"),
    axis.text = element_text(size = 8, color = "#34495E"),
    strip.text = element_text(size = 9, face = "bold", color = "#2C3E50"),
    strip.background = element_rect(fill = "#F8F9FA", color = "#BDC3C7"),
    legend.position = "none",
    panel.grid.major = element_line(color = "#ECF0F1"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# 添加相关系数到散点图
for (i in 1:nrow(cor_results_p3)) {
  group <- cor_results_p3$Group[i]
  corr_val <- round(cor_results_p3$Spearman_corr[i], 3)
  p_val <- ifelse(cor_results_p3$P_value[i] < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_results_p3$P_value[i]))
  
  p_c <- p_c + 
    geom_text(
      data = data.frame(Group = group),
      aes(x = Inf, y = Inf, 
          label = paste("ρ =", corr_val, "\nP =", p_val)),
      hjust = 1.1, vjust = 1.1, size = 2.8, color = "black"
    )
}

# --------------------- Combine All Plots ---------------------
cat("Combining all plots into one figure...\n")

# 创建组合图形布局
combined_plot <- ggdraw() +
  # 上方行：A图（柱状图）和B图（柱状图）
  draw_plot(p_a, 
            x = 0, y = 0.5, 
            width = 0.5, height = 0.5) +
  draw_plot(p_b, 
            x = 0.5, y = 0.5, 
            width = 0.5, height = 0.5) +
  # 下方行：C图（散点图）
  draw_plot(p_c, 
            x = 0, y = 0, 
            width = 1, height = 0.5) +
  # 添加面板标签
  draw_plot_label(
    label = c("A", "B", "C"),
    size = 16,
    fontface = "bold",
    x = c(0, 0.5, 0),
    y = c(1, 1, 0.5)
  )

# 保存组合图形
cat("Saving combined figure...\n")
ggsave("SERPINB1_Combined_Figure.png", 
       combined_plot, 
       width = 16, 
       height = 12, 
       dpi = 300, 
       bg = "white")

ggsave("SERPINB1_Combined_Figure.pdf", 
       combined_plot, 
       width = 16, 
       height = 12, 
       bg = "white")

cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")

cat("Output files created:\n")
cat("1. SERPINB1_Combined_Figure.png - Combined figure with all three plots\n")
cat("2. SERPINB1_Combined_Figure.pdf - High-quality PDF version\n")
cat("3. SERPINB1_Gleason_correlation_results.csv - Correlation results for Part 1\n")
cat("4. SERPINB1_EMAPK_correlation_results.csv - Correlation results for Part 2\n")
cat("5. SERPINB1_Meth_Expr_correlation_results.csv - Correlation results for Part 3\n")