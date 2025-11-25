# ====================== SERPINB1 Methylation-Expression Correlation Analysis ======================
# Objective: Analyze correlation between SERPINB1 promoter methylation and gene expression in cancer/normal tissues
# Pipeline:
#   1. Read methylation data → Extract SERPINB1 promoter CpG sites → Calculate average methylation
#   2. Read mRNA expression data → Extract SERPINB1 expression levels
#   3. Read clinical data → Classify cancer/normal samples
#   4. Integrate all three datasets → Calculate methylation-expression correlation by group
#   5. (Optional) Generate scatter plots by group
# ===========================================================================

# -------------------------- 1. Load Required Packages ----------------------------
library(data.table)    # Efficient large file reading (supports .gz compression)
library(dplyr)         # Data cleaning and pipeline operations
library(tidyr)         # Data reshaping (long/wide format)
library(ggplot2)       # Advanced data visualization
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Provides CpG annotations for 450k array

setwd("Figure2_Data")

# ----------------- 2. Read and Process Methylation Data (450k array) -----------------
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

# Read methylation data (directly decompress .gz file)
methylation_data <- data.table::fread(
  cmd = paste("zcat", meth_file),  # Use zcat to decompress .gz files
  header = TRUE,
  sep = "\t",
  data.table = FALSE,              # Return data.frame
  check.names = FALSE              # Preserve special characters in column names
)

# Data validation and processing
# Structure: Row 1 contains sample IDs (first column empty), Column 1 contains CpG sites (starting from row 2)
if (nrow(methylation_data) < 2) stop("Insufficient rows in methylation data")
if (ncol(methylation_data) < 2) stop("Insufficient columns in methylation data")

# Extract sample IDs (column names excluding first column)
sample_names_meth <- colnames(methylation_data)[-1]

# Extract CpG site names (first column excluding header row)
cpg_names <- methylation_data[-1, 1] %>% as.character()

# Build Beta value matrix (core methylation data)
beta_matrix <- as.matrix(methylation_data[-1, -1])
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names_meth

# Convert to dataframe (preserve CpG site information)
beta_df <- as.data.frame(beta_matrix) %>%
  tibble::rownames_to_column(var = "CpG")

cat("Methylation data loaded:", nrow(beta_matrix), "CpG sites ×", 
    ncol(beta_matrix), "samples\n")

# ----------- 3. Extract SERPINB1 Promoter Region CpG Sites -----------
# Get 450k array annotation (hg19 genome)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter CpG sites in SERPINB1 promoter regions (precise filtering)
cpgs_promoter_serpinb1 <- as.data.frame(anno) %>%
  filter(
    grepl("SERPINB1\\b", UCSC_RefGene_Name),  # Gene annotation contains SERPINB1
    grepl("TSS1500|TSS200|1stExon|5'UTR", UCSC_RefGene_Group)  # Promoter-related regions
  )

cpg_list_promoter <- unique(cpgs_promoter_serpinb1$Name)

# Validate target CpG sites found
if (length(cpg_list_promoter) == 0) {
  stop("Error: No CpG sites found in SERPINB1 promoter regions. Possible reasons:
       1) Gene name misspelled (current: SERPINB1)
       2) Annotation version mismatch (current: hg19)
       3) No promoter probes for this gene on 450k array")
} else {
  cat("Successfully found", length(cpg_list_promoter), "CpG sites in SERPINB1 promoter regions\n")
  cat("Example CpGs:", paste(head(cpg_list_promoter), collapse = ", "), "...\n")
}

# Extract methylation data for target CpGs
meth_promoter <- beta_df %>%
  filter(CpG %in% cpg_list_promoter) %>%
  filter(rowSums(!is.na(select(., -CpG))) > 0)

# Calculate average methylation level per sample
avg_meth_promoter <- meth_promoter %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  group_by(sample_id) %>%
  summarise(
    SERPINB1_Promoter_Beta = mean(Beta_value, na.rm = TRUE),
    CpG_count = n()  # Record number of CpGs used (quality control)
  )

cat("Methylation data processing complete:", nrow(avg_meth_promoter), "samples with promoter average methylation values\n")
cat("CpG count statistics: median=", median(avg_meth_promoter$CpG_count),
    "range=[", min(avg_meth_promoter$CpG_count), "-", 
    max(avg_meth_promoter$CpG_count), "]\n")

# ----------------- 4. Read mRNA Expression Data --------------------

expr_file <- "data/TCGA.PRAD.sampleMap_HiSeqV2.gz"

# Read expression data (preserve original structure)
expr_data <- data.table::fread(
  cmd = paste("zcat", expr_file),  # Decompress .gz file
  header = TRUE,
  sep = "\t",
  data.table = FALSE,
  check.names = FALSE  # Preserve original column names (e.g., sample IDs)
)

# Key modification: Extract gene names (first column, excluding header row)
gene_names <- expr_data[-1, 1] %>% as.character()  # Remove first row, take first column gene names

# Extract sample names (first row, excluding first column empty value/header)
sample_names_expr <- colnames(expr_data)[-1]  # Column names are sample IDs, excluding first column (gene column)

# Build expression matrix (core data: remove first row and first column)
expr_matrix <- expr_data[-1, -1] %>%  # Remove first row (sample ID row) and first column (gene column)
  as.matrix() %>%
  apply(2, as.numeric)  # Ensure all values are numeric (avoid errors from character type)

# Validate row consistency (must match!)
if (length(gene_names) != nrow(expr_matrix)) {
  stop(paste0("Gene name length (", length(gene_names), ") doesn't match expression matrix row count (", nrow(expr_matrix), ")!"))
}

# Assign row names (gene names) and column names (sample IDs)
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- sample_names_expr

cat("Expression data loaded:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "samples\n")

# Convert to dataframe
expr_df <- as.data.frame(expr_matrix) %>%
  tibble::rownames_to_column(var = "Gene")

cat("Expression data loaded:", nrow(expr_matrix), "genes ×", 
    ncol(expr_matrix), "samples\n")

# Extract SERPINB1 expression levels
target_gene <- "SERPINB1"

if (!target_gene %in% expr_df$Gene) {
  stop("Error: Target gene '", target_gene, "' not found. Check gene name spelling")
}

expr_serpinb1 <- expr_df %>%
  filter(Gene == target_gene) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample_id",
    values_to = "SERPINB1_Expression"
  ) %>%
  select(-Gene)  # Gene column no longer needed

cat("Successfully extracted", target_gene, "expression levels for", nrow(expr_serpinb1), "samples\n")

# ----------------- 5. Read and Process Clinical Data --------------------
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

# Read clinical data
clinical_data <- data.table::fread(
  file = clinical_file,  # Note: Clinical data is usually uncompressed
  header = TRUE,
  sep = "\t",
  data.table = FALSE,
  check.names = FALSE
)

# Define key column names (TCGA standard)
sample_id_col <- "sampleID"    # Sample-level ID
sample_type_col <- "sample_type"         # Sample type column

# Validate column existence
if (!sample_id_col %in% colnames(clinical_data)) {
  stop("Error: Sample ID column '", sample_id_col, "' not found. Available columns:",
       paste(colnames(clinical_data), collapse = ", "))
}

if (!sample_type_col %in% colnames(clinical_data)) {
  stop("Error: Sample type column '", sample_type_col, "' not found. Available columns:",
       paste(colnames(clinical_data), collapse = ", "))
}

# Extract key columns and create groups
clinical_df <- clinical_data %>%
  select(
    sample_id = !!sample_id_col,
    sample_type = !!sample_type_col
  ) %>%
  mutate(
    # Create grouping variable: Cancer vs Normal
    Group = case_when(
      grepl("tumor|primary", sample_type, ignore.case = TRUE) ~ "Cancer",
      grepl("normal|solid tissue", sample_type, ignore.case = TRUE) ~ "Normal",
      TRUE ~ NA_character_  # Other types marked as NA
    )
  ) %>%
  filter(!is.na(Group))  # Keep only clearly grouped samples

# Group statistics
group_counts <- clinical_df %>%
  dplyr::count(Group)
cat("\nClinical data group statistics:\n")
print(group_counts)

# ----------------- 6. Data Integration and Merging --------------------
# Three-step merge: Methylation + Expression + Grouping
merged_data <- avg_meth_promoter %>%
  inner_join(expr_serpinb1, by = "sample_id") %>%
  inner_join(clinical_df, by = "sample_id")

# Validate merge results
if (nrow(merged_data) == 0) {
  stop("Error: No valid samples after data merging. Check sample ID matching")
}

cat("\nData integration complete:", nrow(merged_data), "samples available for analysis\n")
cat("Final group distribution:\n")
print(table(merged_data$Group))

write.csv(merged_data, "SERPINB1_methylation_vs_mRNA_merged_data.csv")

# ----------------- 7. Correlation Analysis --------------------
# Calculate methylation-expression correlation by group
cor_results <- merged_data %>%
  group_by(Group) %>%
  summarise(
    Spearman_corr = cor(SERPINB1_Promoter_Beta, SERPINB1_Expression, 
                        method = "spearman", use = "complete.obs"),
    P_value = suppressWarnings(  # Suppress small sample warnings
      cor.test(SERPINB1_Promoter_Beta, SERPINB1_Expression, 
               method = "spearman", exact = FALSE)$p.value
    ),
    N = n(),
    .groups = "drop"
  )

# Results output
cat("\n===== SERPINB1 Promoter Methylation-Expression Correlation Results =====\n")
print(cor_results)
write.csv(cor_results, "correlation between SERPINB1 promoter methylation and gene expression in cancer_normal tissues.csv")
# ----------------- 8. Visualization (Optional) --------------------
# Create scatter plot
p <- ggplot(merged_data, 
            aes(x = SERPINB1_Promoter_Beta, 
                y = SERPINB1_Expression, 
                color = Group)) +
  geom_point(alpha = 0.6, size = 2.5) +  # Semi-transparent points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +  # Add trend line
  facet_wrap(~ Group, scales = "free") +  # Show by group
  labs(
    title = "SERPINB1 Promoter Methylation-Gene Expression Correlation",
    subtitle = paste("TCGA-PRAD Data | Cancer samples:", 
                     sum(merged_data$Group == "Cancer"),
                     "| Normal samples:", 
                     sum(merged_data$Group == "Normal")),
    x = "Promoter Region Average Methylation (β-value)",
    y = paste(target_gene, "Gene Expression"),
    caption = paste("Analysis includes", length(cpg_list_promoter), "promoter CpG sites")
  ) +
  scale_color_manual(
    values = c("Cancer" = "#E15759", "Normal" = "#4E79A7")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",  # Remove legend (facets differentiate)
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

# Add correlation coefficient annotations
for (i in 1:nrow(cor_results)) {
  group <- cor_results$Group[i]
  corr_val <- round(cor_results$Spearman_corr[i], 3)
  p_val <- ifelse(cor_results$P_value[i] < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_results$P_value[i]))
  
  p <- p + 
    geom_text(
      data = data.frame(Group = group),
      aes(x = Inf, y = Inf, 
          label = paste("ρ =", corr_val, "\nP =", p_val)),
      hjust = 1.1, vjust = 1.1, size = 4.5, color = "black"
    )
}

# Save high-quality image
ggsave(paste0("SERPINB1_Methylation_vs_Expression_", Sys.Date(), ".pdf"),
       plot = p, width = 10, height = 5, dpi = 300)

# Print completion message
cat("\n===== Analysis Complete =====\n")
cat("Results saved as PDF file\n")

