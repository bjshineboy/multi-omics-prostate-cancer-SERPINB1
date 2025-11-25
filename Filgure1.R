# ===================== Methylation Analysis Pipeline =====================
# Objective: Analyze methylation differences of SERPINB1-related CpG sites between normal and tumor tissues in prostate cancer
# Data Source: TCGA-PRAD Methylation 450K array data
# Author: Bioinformatics Analysis Team
# Date: 2025-06-27
# ========================================================================

# --------------------- 1. Load Required Packages ---------------------
library(data.table)    # Efficient reading of large data files
library(dplyr)         # Data manipulation and piping
library(tidyr)         # Data transformation (wide to long)
library(ggplot2)       # Advanced data visualization
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # 450K array hg19 annotation
library(stringr)       # String processing
library(tibble)        # Provides rownames_to_column function
library(patchwork)     # Plot composition

setwd("Figure1_Data")
# --------------------- 2. Read Methylation Data -----------------
# File paths
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

# Read compressed methylation data efficiently
cat("Reading methylation data...\n")
methylation_data <- fread(
  cmd = paste("zcat", meth_file),
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

# Data verification and preprocessing
cat("Methylation data dimensions:", dim(methylation_data), "\n")
cat("Methylation values of first 5 CpG sites in first 5 samples:\n")
print(head(methylation_data[, 1:5]))

# 2. Extract and set row and column names
# Sample names
sample_names <- colnames(methylation_data)[-1]  # Remove first column

# CpG site names (first column without header)
cpg_names <- as.character(methylation_data[-1, 1])

# Beta matrix (remove first row and column)
beta_matrix <- as.matrix(methylation_data[-1, -1])

# Verify dimension matching
if (length(cpg_names) != nrow(beta_matrix)) {
  stop(paste("Dimension mismatch: cpg_names length =", length(cpg_names), 
             ", beta_matrix rows =", nrow(beta_matrix)))
}

# Set row and column names
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names

# --------------------- 3. Get SERPINB1-related CpG Sites ---------

# Get SERPINB1-related CpG site annotations (corrected column names)
cat("\nGetting annotations for SERPINB1 gene-related CpG sites...\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter SERPINB1-related CpG sites (exact match gene symbol)
cpgs_serpinb1 <- as.data.frame(anno) %>%
  dplyr::filter(str_detect(UCSC_RefGene_Name, "SERPINB1\\b"))  # Use word boundary for exact match

# Extract detailed CpG site information (using correct column names)
cpgs_serpinb1_info <- cpgs_serpinb1 %>%
  dplyr::select(
    Name,           # CpG site name (keep)
    chr,            # Chromosome (corrected: original CHR → chr)
    pos,            # Genomic position (corrected: original MAPINFO → pos)
    UCSC_RefGene_Name,  # Associated gene (keep)
    UCSC_RefGene_Group  # Gene functional region (keep)
  ) %>%
  arrange(chr, pos)  # Sort by chromosome and position

# Verify results
cat("Found", nrow(cpgs_serpinb1_info), "CpG sites directly related to SERPINB1\n")
if(nrow(cpgs_serpinb1_info) > 0) {
  cat("Detailed information for first 5 SERPINB1-related CpG sites:\n")
  print(head(cpgs_serpinb1_info, 5))
  
  cpg_list <- cpgs_serpinb1_info$Name
} else {
  stop("No SERPINB1-related CpG sites found, please check gene name or annotation database")
}

# --------------------- 4. Extract Target CpG Site Data -------------
# Filter methylation data for target CpG sites
target_beta <- beta_matrix[rownames(beta_matrix) %in% cpg_list, ]

# Convert to long format dataframe
beta_long <- as.data.frame(target_beta) %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  mutate(Beta_value = as.numeric(Beta_value))  # Ensure numeric type

cat("\nTarget CpG site data summary:\n")
cat("Number of samples:", length(unique(beta_long$sample_id)), "\n")
cat("Number of CpG sites:", length(unique(beta_long$CpG)), "\n")
cat("Total observations:", nrow(beta_long), "\n")

# --------------------- 5. Read Clinical Data -------------------
cat("\nReading clinical data...\n")
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

clinical_data <- fread(
  clinical_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE,
  na.strings = c("", "NA", "[Not Available]")
)

# Clinical data verification
cat("Clinical data dimensions:", dim(clinical_data), "\n")
cat("First 5 clinical data columns:\n")
print(head(colnames(clinical_data), 5))

# Identify sample ID column (common in TCGA data)
id_cols <- c("sampleID", "bcr_patient_barcode", "bcr_sample_barcode")
sample_id_col <- intersect(id_cols, colnames(clinical_data))[1]

if(is.na(sample_id_col)) {
  stop("Sample ID column not found, please check clinical data column names")
}
cat("Using column '", sample_id_col, "' as sample ID\n")

# Extract key clinical information
clinical_df <- clinical_data %>%
  dplyr::select(
    sample_id = !!sample_id_col,
    sample_type = "sample_type"
  ) %>%
  mutate(
    # Create standardized grouping variable - modified group names to Tumor and Normal
    group = case_when(
      grepl("Primary Tumor", sample_type, ignore.case = TRUE) ~ "Tumor",
      grepl("Solid Tissue Normal", sample_type, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(group %in% c("Tumor", "Normal"))  # Keep only key groups

# Group sample statistics
cat("\nClinical group sample counts:\n")
print(table(clinical_df$group))
# Save merged data
write.csv(clinical_df, "clinical_df_data.csv", row.names = FALSE)
# --------------------- 6. Merge Methylation and Clinical Data -----------
meth_clin <- beta_long %>%
  inner_join(clinical_df, by = "sample_id") %>%
  left_join(cpgs_serpinb1_info, by = c("CpG" = "Name"))  # Add annotation information

# Add simplified labels for visualization
meth_clin <- meth_clin %>%
  mutate(
    Position = paste0(chr, ":", pos),  # Use lowercase chr and pos
    Gene_Region = case_when(
      grepl("TSS", UCSC_RefGene_Group) ~ "Promoter",
      grepl("Body", UCSC_RefGene_Group) ~ "Gene Body",
      grepl("5'UTR", UCSC_RefGene_Group) ~ "5'UTR",
      grepl("3'UTR", UCSC_RefGene_Group) ~ "3'UTR",
      TRUE ~ "Other Region"
    )
  )

cat("\nMerged data summary:\n")
print(dim(meth_clin))
cat("Group sample counts:\n")
print(table(meth_clin$group))

# Simplify position labels (avoid long CpG names)
position_labels <- meth_clin %>%
  distinct(CpG, Position) %>%
  mutate(Label = paste0("CpG-", seq_along(CpG), " (", Position, ")"))

meth_clin <- meth_clin %>%
  dplyr::filter(!is.na(Beta_value)) %>%
  left_join(position_labels, by = c("CpG", "Position"))

# Save merged data
write.csv(meth_clin, "merged_methylation_clinical_data.csv", row.names = FALSE)
# --------------------- 7. Data Visualization ---------------------
cat("\nGenerating methylation level visualizations...\n")

# 7.1 Boxplot grouped by CpG site (modified group labels to Tumor and Normal)
p1 <- ggplot(meth_clin, aes(x = group, y = Beta_value, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ CpG, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Tumor" = "#4E79A7", "Normal" = "#F28E2B")) +
  labs(
    title = "Methylation Levels of SERPINB1-related CpG Sites",
    subtitle = "Prostate Cancer Tumor Tissue vs Normal Tissue",
    x = "Sample Type",
    y = "Methylation β-value",
    caption = paste("Data Source: TCGA-PRAD | Number of CpG sites:", length(unique(meth_clin$CpG)))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p1)

# 7.2 Heatmap grouped by genomic position (modified group labels to Tumor and Normal)


p2 <- ggplot(meth_clin, aes(x = sample_id, y = CpG, fill = Beta_value)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(
    low = "#4E79A7", 
    mid = "white", 
    high = "#F28E2B",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  facet_grid(Gene_Region ~ group, scales = "free", space = "free") +
  labs(
    title = "Heatmap of SERPINB1-related CpG Site Methylation",
    subtitle = "Grouped by Genomic Region and Sample Type",
    x = "Sample",
    y = "CpG Site (Chromosomal Position)",
    fill = "Methylation β-value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p2)

# 7.3 Methylation levels grouped by genomic region (modified group labels to Tumor and Normal)
p3 <- ggplot(meth_clin, aes(x = Gene_Region, y = Beta_value, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 23, 
    size = 3, 
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = c("Tumor" = "#4E79A7", "Normal" = "#F28E2B")) +
  labs(
    title = "Methylation Levels by Genomic Region",
    x = "Genomic Region",
    y = "Methylation β-value",
    fill = "Sample Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p3)

# --------------------- 8. Statistical Analysis -----------------------
cat("\nPerforming statistical analysis...\n")

# 8.1 Group t-test by CpG site (modified group labels to Tumor and Normal)
results <- meth_clin %>%
  group_by(CpG, Position, Gene_Region) %>%
  summarize(
    Tumor_Mean = mean(Beta_value[group == "Tumor"], na.rm = TRUE),
    Normal_Mean = mean(Beta_value[group == "Normal"], na.rm = TRUE),
    Mean_Difference = Tumor_Mean - Normal_Mean,
    p_value = tryCatch(
      t.test(Beta_value ~ group)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Output significant results
sig_results <- results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

cat("\nSignificantly different CpG sites:\n")
if(nrow(sig_results) > 0) {
  print(sig_results)
  write.csv(sig_results, "SERPINB1_methylation_cpg_sig_results.csv", row.names = FALSE)
  
} else {
  cat("No significantly different CpG sites found (p < 0.05)\n")
}

# 8.2 Overall difference analysis (all CpG sites combined)
overall_test <- t.test(Beta_value ~ group, data = meth_clin)
cat("\nOverall methylation level difference (all CpG sites):\n")
print(overall_test)
write.csv( overall_test$estimate, "SERPINB1_methylation_overall_test_estimate.csv", row.names = FALSE)
# --------------------- 9. Save Results -----------------------
cat("\nSaving analysis results...\n")

# Save statistical results
write.csv(results, "SERPINB1_methylation_statistical_results.csv", row.names = FALSE)

# Combine plots
combined_plot <- p1 + p3 + p2 + 
  plot_layout(
    nrow = 2,
    heights = c(1, 1.5),
    design = "
    AB
    CC
    "
  ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.')



# Save visualizations
ggsave("SERPINB1_methylation_boxplot.png", p1, width = 12, height = 8, dpi = 300)
ggsave("SERPINB1_methylation_heatmap.png", p2, width = 14, height = 10, dpi = 300)
ggsave("SERPINB1_region_methylation.png", p3, width = 10, height = 6, dpi = 300)
ggsave("SERPINB1_combined_plot_methylation.png", combined_plot, width = 16, height = 12, dpi = 300)
cat("\nAnalysis complete! Results saved to current working directory\n")


##########
library(dplyr)
library(tidyr)

summary_mean <- meth_clin %>%
  group_by(Gene_Region, group) %>%
  summarise(
    Beta_mean = mean(Beta_value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(Beta_mean, n),
    names_sep = "_"
  )

 
summary_mean <- summary_mean %>%
  mutate(
    Beta_diff = Beta_mean_Tumor - Beta_mean_Normal
  )

# 
p_values <- meth_clin %>%
  filter(!is.na(group) & !is.na(Beta_value)) %>%
  group_by(Gene_Region) %>%
  summarise(
    p_value = t.test(Beta_value ~ group, data = .)$p.value,
    .groups = 'drop'
  )

# 
final_stats <- summary_mean %>%
  left_join(p_values, by = "Gene_Region")

# 
final_stats <- final_stats %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


write.csv(final_stats, "SERPINB1_methylation_gene_region.csv", row.names = FALSE)


