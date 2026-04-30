# ===================== Methylation Analysis Pipeline =====================
# Objective: Analyze methylation differences of SERPINB1-related CpG sites
# Data Source: TCGA-PRAD Methylation 450K array data
# Author: Bioinformatics Analysis Team
# Date: 2025-06-27
# ========================================================================

# --------------------- 1. Load Required Packages ---------------------
library(data.table)    # Efficient reading of large data files
library(dplyr)         # Data manipulation and piping
library(tidyr)         # Data transformation
library(ggplot2)       # Advanced data visualization
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # 450K array hg19 annotation
library(stringr)       # String processing
library(tibble)        # Provides rownames_to_column function
library(patchwork)     # Plot composition
library(ggpubr)        # Publication-ready plots
library(RColorBrewer)  # Color palettes
library(forcats)       # Factor manipulation
library(scales)        # Scale functions
library(ggrepel)       # For better text labels
library(viridis)       # Colorblind-friendly color schemes
library(ComplexHeatmap) # For better heatmaps
library(circlize)      # Color scales for heatmaps
library(ggsci)         # Scientific journal color palettes
library(gridExtra)     # For table generation
library(grid)          # For low-level grid functions

setwd("/work/singleCell/prostate_singecell_20250627/GSE_Raw/sci/SERPNB1_sci_20250829/code/2/Figure1_Data")

# Define color palette
my_colors <- c("Normal" = "#2E8B57",  # Sea Green
               "Tumor" = "#DC143C")   # Crimson Red

# --------------------- 2. Read Methylation Data -----------------
# File paths
meth_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

cat("Reading methylation data...\n")
methylation_data <- fread(
  cmd = paste("zcat", meth_file),
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

cat("Methylation data dimensions:", dim(methylation_data), "\n")

# Extract and set row and column names
sample_names <- colnames(methylation_data)[-1]
cpg_names <- methylation_data[[1]][-1]
beta_matrix <- as.matrix(methylation_data[-1, -1])

rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names

cat("Beta matrix dimensions:", dim(beta_matrix), "\n")

# --------------------- 3. Get SERPINB1-related CpG Sites ---------
cat("\nGetting annotations for SERPINB1 gene-related CpG sites...\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter SERPINB1-related CpG sites
cpgs_serpinb1 <- as.data.frame(anno) %>%
  dplyr::filter(str_detect(UCSC_RefGene_Name, "SERPINB1\\b"))

cpgs_serpinb1_info <- cpgs_serpinb1 %>%
  dplyr::select(
    Name,
    chr,
    pos,
    UCSC_RefGene_Name,
    UCSC_RefGene_Group
  ) %>%
  arrange(chr, pos)

cat("Found", nrow(cpgs_serpinb1_info), "CpG sites related to SERPINB1\n")

if(nrow(cpgs_serpinb1_info) > 0) {
  cpg_list <- cpgs_serpinb1_info$Name
} else {
  stop("No SERPINB1-related CpG sites found")
}

# --------------------- 4. Extract Target CpG Site Data -------------
# Filter methylation data for target CpG sites
target_rows <- rownames(beta_matrix) %in% cpg_list
cat("Found", sum(target_rows), "matching CpG sites in data\n")

if(sum(target_rows) == 0) {
  stop("No matching CpG sites found in methylation data")
}

target_beta <- beta_matrix[target_rows, ]

# Convert to long format
beta_long <- as.data.frame(target_beta) %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  mutate(Beta_value = as.numeric(Beta_value))

cat("\nTarget CpG data summary:\n")
cat("Number of samples:", length(unique(beta_long$sample_id)), "\n")
cat("Number of CpG sites:", length(unique(beta_long$CpG)), "\n")

# --------------------- 5. Read Clinical Data -------------------
cat("\nReading clinical data...\n")
clinical_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

clinical_data <- fread(
  clinical_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE,
  na.strings = c("", "NA", "[Not Available]", "[Not Applicable]")
)

cat("Clinical data dimensions:", dim(clinical_data), "\n")

# Create clinical dataframe
clinical_df <- clinical_data %>%
  dplyr::select(
    sample_id = sampleID,
    sample_type = sample_type
  ) %>%
  mutate(
    group = case_when(
      grepl("Primary Tumor|Tumor|Primary", sample_type, ignore.case = TRUE) ~ "Tumor",
      grepl("Solid Tissue Normal|Normal|Control", sample_type, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(group %in% c("Tumor", "Normal"))

cat("\nClinical group distribution:\n")
print(table(clinical_df$group))

# --------------------- 6. Merge Methylation and Clinical Data -----------
# Clean sample IDs for matching
clean_sample_id <- function(sample_id) {
  substr(sample_id, 1, 15)
}

beta_long <- beta_long %>%
  mutate(sample_clean = clean_sample_id(sample_id))

clinical_df <- clinical_df %>%
  mutate(sample_clean = clean_sample_id(sample_id))

# Merge data
meth_clin <- beta_long %>%
  inner_join(clinical_df, by = "sample_clean") %>%
  left_join(cpgs_serpinb1_info, by = c("CpG" = "Name"))

cat("\nMerged data dimensions:", dim(meth_clin), "\n")
cat("Group distribution in merged data:\n")
print(table(meth_clin$group))

# Filter out NA values
meth_clin <- meth_clin %>%
  dplyr::filter(!is.na(Beta_value))

cat("Data after removing NA:", nrow(meth_clin), "rows\n")

# Create simplified labels
meth_clin <- meth_clin %>%
  mutate(
    Position = paste0(chr, ":", pos),
    Gene_Region = case_when(
      grepl("TSS", UCSC_RefGene_Group) ~ "Promoter",
      grepl("Body", UCSC_RefGene_Group) ~ "Gene Body",
      grepl("5'UTR", UCSC_RefGene_Group) ~ "5'UTR",
      grepl("3'UTR", UCSC_RefGene_Group) ~ "3'UTR",
      grepl("1stExon", UCSC_RefGene_Group) ~ "First Exon",
      TRUE ~ "Other"
    ),
    Gene_Region = factor(Gene_Region, 
                         levels = c("Promoter", "5'UTR", "First Exon", 
                                    "Gene Body", "3'UTR", "Other"))
  )

# Create CpG labels with direct CpG IDs
cpg_labels <- meth_clin %>%
  distinct(CpG, chr, pos, Gene_Region) %>%
  arrange(chr, pos) %>%
  mutate(
    CpG_label = CpG,  # Use direct CpG ID
    Position_short = paste0(chr, ":", pos),
    CpG_label = factor(CpG_label, levels = unique(CpG_label))  # Keep as factor for ordering
  )

meth_clin <- meth_clin %>%
  left_join(cpg_labels, by = c("CpG", "chr", "pos", "Gene_Region"))

# Save data
write.csv(meth_clin, "merged_methylation_clinical_data.csv", row.names = FALSE)
cat("\nSaved merged data to: merged_methylation_clinical_data.csv\n")

# --------------------- 7. Statistical Analysis -----------------------
cat("\nPerforming statistical analysis...\n")

# 7.1 T-test by CpG site
results <- meth_clin %>%
  group_by(CpG, Position, Gene_Region, CpG_label) %>%
  summarise(
    Tumor_Mean = mean(Beta_value[group == "Tumor"], na.rm = TRUE),
    Normal_Mean = mean(Beta_value[group == "Normal"], na.rm = TRUE),
    Tumor_SD = sd(Beta_value[group == "Tumor"], na.rm = TRUE),
    Normal_SD = sd(Beta_value[group == "Normal"], na.rm = TRUE),
    Mean_Difference = Tumor_Mean - Normal_Mean,
    p_value = tryCatch(
      t.test(Beta_value ~ group)$p.value,
      error = function(e) NA
    ),
    n_Tumor = sum(group == "Tumor"),
    n_Normal = sum(group == "Normal"),
    .groups = "drop"
  ) %>%
  mutate(
    p_adjust = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      p_adjust < 0.001 ~ "***",
      p_adjust < 0.01 ~ "**",
      p_adjust < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adjust)

# Save statistical results
write.csv(results, "SERPINB1_methylation_statistical_results.csv", row.names = FALSE)

# --------------------- 8. Data Visualization ---------------------
cat("\nGenerating visualizations...\n")

# Define a common theme with optimized font sizes
common_theme <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

# 8.1 Create Figure A: Boxplot for Top CpG Sites
cat("\nCreating Figure A: Boxplot for top CpG sites...\n")

# Select top 10 CpG sites for display (based on significance or mean difference)
top_cpgs <- results %>%
  arrange(desc(abs(Mean_Difference))) %>%
  head(10) %>%
  pull(CpG)

# Filter data for top CpG sites
top_cpg_data <- meth_clin %>%
  filter(CpG %in% top_cpgs) %>%
  mutate(
    CpG_label = factor(CpG, levels = top_cpgs)
  )

# Create boxplot for top CpG sites
p_boxplot_cpg <- ggplot(top_cpg_data, aes(x = group, y = Beta_value, fill = group)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.85,
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.color = "gray30",
    outlier.size = 1,
    lwd = 0.3
  ) +
  scale_fill_manual(values = my_colors, name = "Sample Type") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  facet_wrap(~ CpG_label, ncol = 5, scales = "free_y") +
  labs(
    title = "A: Differential Methylation of Top SERPINB1 CpG Sites",
    x = "Sample Type",
    y = "DNA Methylation (β-value)"
  ) +
  common_theme(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8),
    plot.title = element_text(size = 10),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top"
  )

# Add statistical comparisons to each CpG
p_boxplot_cpg <- p_boxplot_cpg + 
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    label.y = 0.95,
    size = 3,
    bracket.size = 0.3,
    tip.length = 0.01
  )

# 8.2 Create Figure B: Volcano Plot
cat("\nCreating Figure B: Volcano plot...\n")
# 8.2 Create Figure B: Volcano Plot
cat("\nCreating Figure B: Volcano plot...\n")

volcano_data <- results %>%
  mutate(
    log10_p = -log10(p_adjust),
    is_significant = p_adjust < 0.05
  )

p_volcano <- ggplot(volcano_data, aes(x = Mean_Difference, y = log10_p, 
                                      color = is_significant, 
                                      alpha = is_significant)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#D55E00"), guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1), guide = "none") +
  # 修改这里：标注所有显著的CpG位点，而不仅仅是前5个
  geom_text_repel(
    data = volcano_data %>% filter(p_adjust < 0.05),  # 移除了head(5)限制
    aes(label = CpG),
    size = 2.5,
    max.overlaps = 50,  # 增加最大重叠数
    min.segment.length = 0.1,  # 减小最小段长度
    box.padding = 0.5,  # 增加标签框的内边距
    point.padding = 0.3,  # 增加点与标签之间的距离
    segment.size = 0.2,  # 减小连接线的大小
    segment.color = "grey50",  # 设置连接线颜色
    segment.alpha = 0.5,  # 设置连接线透明度
    direction = "both",  # 允许标签在两个方向上移动
    force = 1,  # 增加排斥力
    nudge_x = 0.1,  # 水平方向微调
    nudge_y = 0.1,  # 垂直方向微调
    max.time = 1,  # 设置最大计算时间
    max.iter = 10000  # 增加最大迭代次数
  ) +
  # 为不显著的位点添加更小的、透明的标签（可选）
  geom_text_repel(
    data = volcano_data %>% filter(p_adjust >= 0.05),
    aes(label = CpG),
    size = 1.5,  # 更小的字体
    max.overlaps = 20,
    min.segment.length = 0.2,
    alpha = 0.3,  # 透明度
    color = "gray40",  # 灰色标签
    segment.size = 0.1,
    segment.alpha = 0.2
  ) +
  labs(
    title = "B: Volcano Plot: Differential Methylation of SERPINB1 CpG Sites",
    x = "Methylation Difference (β-Tumor - β-Normal)",
    y = expression(-log[10]("FDR-adjusted p-value"))
  ) +
  common_theme(base_size = 9) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10)
  )
# 8.3 Create Figure C: Heatmap by Genomic Region
cat("\nCreating Figure C: Heatmap by genomic region...\n")

# Prepare data for heatmap
heatmap_data <- meth_clin %>%
  group_by(CpG_label, group, Gene_Region) %>%
  summarise(
    mean_beta = mean(Beta_value, na.rm = TRUE),
    .groups = "drop"
  )

p_heatmap <- ggplot(heatmap_data, aes(x = group, y = CpG_label, fill = mean_beta)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_beta)), 
            size = 2.5, color = "black") +
  scale_fill_viridis_c(
    option = "plasma",
    direction = -1,
    limits = c(0, 1),
    name = "Mean Methylation (β)"
  ) +
  facet_grid(Gene_Region ~ ., scales = "free_y", space = "free", switch = "y") +
  labs(
    title = "C: SERPINB1 Methylation Profile by Genomic Region",
    x = "Sample Type",
    y = "CpG Sites"
  ) +
  common_theme(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 8),
    legend.position = "right",
    legend.key.height = unit(1, "cm"),
    plot.title = element_text(size = 10)
  )

# 8.4 Create Figure D: Boxplot by Genomic Region
cat("\nCreating Figure D: Boxplot by genomic region...\n")

p_boxplot_region <- ggplot(meth_clin, aes(x = Gene_Region, y = Beta_value, fill = group)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.85,
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.color = "gray30",
    outlier.size = 1,
    lwd = 0.3
  ) +
  scale_fill_manual(values = my_colors, name = "Sample Type") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "D: Methylation Distribution by Genomic Region",
    x = "Genomic Region",
    y = "DNA Methylation (β-value)"
  ) +
  common_theme(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(size = 10)
  )

# Add statistical comparisons
p_boxplot_region <- p_boxplot_region + 
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    label.y = 0.95,
    size = 3,
    bracket.size = 0.3,
    tip.length = 0.01
  )

# --------------------- 9. Create Final Combined Plot -----------
cat("\nCreating combined figure...\n")

# Convert ggplot objects to grobs
boxplot_cpg_grob <- ggplotGrob(p_boxplot_cpg)
volcano_grob <- ggplotGrob(p_volcano)
heatmap_grob <- ggplotGrob(p_heatmap)
boxplot_region_grob <- ggplotGrob(p_boxplot_region)

# Create a 2x2 layout
combined_plot <- grid.arrange(
  # First row: Boxplot for top CpG sites (Figure A) and Volcano plot (Figure B)
  arrangeGrob(
    boxplot_cpg_grob,
    volcano_grob,
    ncol = 2,
    widths = c(1.2, 0.8),  # Adjusted to give more space to boxplot
    padding = unit(5, "mm")
  ),
  # Second row: Heatmap (Figure C) and Boxplot by region (Figure D)
  arrangeGrob(
    heatmap_grob,
    boxplot_region_grob,
    ncol = 2,
    widths = c(1, 1),
    padding = unit(5, "mm")
  ),
  nrow = 2,
  heights = c(1, 1.2),  # Adjusted to give more space to bottom row
  # Add main title
  top = textGrob(
    "SERPINB1 Methylation Analysis in Prostate Cancer",
    gp = gpar(fontsize = 16, fontface = "bold"),
    just = "center"
  ),
  # Add footer
  bottom = textGrob(
    paste(
      "Data source: TCGA-PRAD |",
      "Total CpG sites:", nrow(results), "|",
      "Tumor samples:", sum(clinical_df$group == "Tumor"), "|",
      "Normal samples:", sum(clinical_df$group == "Normal")
    ),
    gp = gpar(fontsize = 9, col = "gray40"),
    just = "center"
  )
)

# Save combined plot
cat("\nSaving combined figure...\n")

# Save as PDF (for publication quality)
ggsave(
  "Figure5_SERPINB1_combined_analysis.pdf",
  combined_plot,
  width = 12,  # Slightly wider to accommodate the boxplot
  height = 8,  # Slightly taller to accommodate the boxplot
  device = "pdf"
)

# Save as PNG
ggsave(
  "Figure5_SERPINB1_combined_analysis.png",
  combined_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

 
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
 
cat("\nGenerated output files:\n")
cat("1. merged_methylation_clinical_data.csv - Merged methylation and clinical data\n")
cat("2. SERPINB1_methylation_statistical_results.csv - Statistical results\n")
cat("3. Figure5_SERPINB1_combined_analysis.pdf - Combined figure (PDF)\n")
cat("4. Figure5_SERPINB1_combined_analysis.png - Combined figure (PNG)\n")
cat("\nFigure layout:\n")
cat("- Figure A (Top-left): Boxplot of top differentially methylated CpG sites\n")
cat("- Figure B (Top-right): Volcano plot of differential methylation\n")
cat("- Figure C (Bottom-left): Heatmap of methylation by genomic region\n")
cat("- Figure D (Bottom-right): Boxplot of methylation distribution by genomic region\n")

