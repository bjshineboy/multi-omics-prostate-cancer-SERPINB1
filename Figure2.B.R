# ==============================================================================
# Prostate Cancer Methylation Data Analysis: Correlation between SERPINB1-related CpG sites and Gleason Score
# 
# Objective: Analyze the correlation between methylation levels of SERPINB1-related CpG sites
#            and Gleason scores in TCGA-PRAD data
#
# Pipeline Overview:
#   1. Load required R packages
#   2. Read and preprocess methylation data
#   3. Extract Gleason score clinical data
#   4. Identify SERPINB1-related CpG sites
#   5. Integrate methylation data with clinical data
#   6. Calculate Spearman correlations
#   7. Visualize results (heatmap)
#
# Input Files:
#   - Methylation data: TCGA.PRAD.sampleMap_HumanMethylation450.gz
#   - Clinical data: TCGA.PRAD.sampleMap_PRAD_clinicalMatrix
#
# Output:
#   - Correlation statistics table
#   - Heatmap visualization
#
# Author: Bioinformatics Analysis Team
# Date: 2025-06-27
# ==============================================================================

# --------------------- 1. Load Required Packages ---------------------
library(data.table)     # Efficient reading of large data files
library(dplyr)          # Data manipulation and piping
library(tidyr)          # Data transformation (wide to long)
library(ggplot2)        # Basic visualization
library(ComplexHeatmap)  # Advanced heatmap plotting
library(circlize)       # Heatmap color configuration
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # 450K array annotation
library(dplyr)
library(tidyr)
library(tibble) 

setwd("Figure2_Data")
# --------------------- 2. Read and Preprocess Methylation Data ---------------------
# Methylation data file path
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

cat("Step 1: Reading methylation data...\n")
# Use fread to efficiently read gzipped file
methylation_data <- tryCatch(
  {
    fread(
      cmd = paste("zcat", meth_file),
      header = TRUE,
      sep = "\t",
      data.table = FALSE,
      showProgress = TRUE  # Show reading progress
    )
  },
  error = function(e) {
    stop(paste("Failed to read methylation data:", e$message))
  }
)

# Data validation
if (nrow(methylation_data) < 2 || ncol(methylation_data) < 2) {
  stop("Incorrect methylation data format: rows=", nrow(methylation_data), ", columns=", ncol(methylation_data))
}

cat("Methylation data dimensions:", dim(methylation_data), "\n")
cat("Data preview: first 5 rows and 5 columns\n")
print(head(methylation_data[, 1:5], 3))

# Extract sample names (row 1, columns 2+)
sample_names <- colnames(methylation_data)[-1]
cat("Number of samples:", length(sample_names), "\n")

# Extract CpG site names (column 1, rows 2+)
cpg_names <- as.character(methylation_data[-1, 1])

# Extract Beta value matrix and convert to numeric matrix
beta_matrix <- as.matrix(methylation_data[-1, -1])
mode(beta_matrix) <- "numeric"  # Ensure numeric type

# Set row and column names
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names

# Convert to dataframe format (keep CpG sites as column)
beta_df <- as.data.frame(beta_matrix) %>%
  rownames_to_column(var = "CpG")

cat("Methylation matrix preprocessing complete: CpG sites=", nrow(beta_matrix), ", samples=", ncol(beta_matrix), "\n")

# Save intermediate methylation data
#write.csv(beta_df, "methylation_beta_values_all_cpgs.csv", row.names = FALSE)
cat("Saved methylation beta values to: methylation_beta_values_all_cpgs.csv\n")

# --------------------- 3. Read Clinical Data and Extract Gleason Scores ---------------------
# Clinical data file path
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

cat("Step 2: Reading clinical data and extracting Gleason scores...\n")
# Read clinical data
clinical_data <- tryCatch(
  {
    fread(
      clinical_file,
      header = TRUE,
      sep = "\t",
      data.table = FALSE,
      na.strings = c("", "NA", "[Not Available]")  # Specify missing value representations
    )
  },
  error = function(e) {
    stop(paste("Failed to read clinical data:", e$message))
  }
)

# Data validation
if (nrow(clinical_data) < 1) {
  stop("Clinical data is empty")
}

cat("Clinical data dimensions:", dim(clinical_data), "\n")
cat("Clinical data column names:", head(colnames(clinical_data)), "...\n")


clinical_sample_id_col <- "sampleID"
gleason_col <- "gleason_score"  # Use first match

# Extract sample IDs and Gleason scores
clinical_df <- clinical_data %>%
  select(
    sample_id = !!clinical_sample_id_col,
    Gleason = !!gleason_col,
    sample_type
  ) %>%
  mutate(
    # Convert Gleason score to numeric
    Gleason = as.numeric(as.character(Gleason)),
  ) %>%
  filter(
    !is.na(Gleason),          # Remove missing values
    sample_type == "Primary Tumor",
    Gleason >= 5, Gleason <= 10  # Validate valid range
  )

# Validate Gleason score distribution
cat("Gleason score distribution:\n")
print(table(clinical_df$Gleason))

# Sample type statistics
cat("Sample type distribution:\n")
print(table(clinical_df$sample_type))

cat("Extracted Gleason scores for", nrow(clinical_df), "samples\n")

# Save clinical data
write.csv(clinical_df, "clinical_data_with_gleason_scores.csv", row.names = FALSE)
cat("Saved clinical data to: clinical_data_with_gleason_scores.csv\n")

# --------------------- 4. Identify SERPINB1-related CpG Sites ---------------------
cat("Step 3: Identifying SERPINB1-related CpG sites...\n")

# Get methylation array annotation data
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter SERPINB1-related CpG sites
cpgs_serpinb1 <- tryCatch(
  {
    as.data.frame(anno) %>%
      filter(grepl("SERPINB1\\b", UCSC_RefGene_Name))  # Exact gene name match
  },
  error = function(e) {
    stop("Failed to filter SERPINB1-related CpG sites:", e$message)
  }
)

# Check filtering results
if (nrow(cpgs_serpinb1) == 0) {
  stop("No CpG sites related to SERPINB1 found, please check gene name or annotation database")
}

cpg_list <- cpgs_serpinb1$Name
cat("Found", length(cpg_list), "CpG sites related to SERPINB1\n")

# Extract CpG site details
cpgs_serpinb1_info <- cpgs_serpinb1 %>%
  select(
    Name, chr, pos, 
    UCSC_RefGene_Name, 
    UCSC_RefGene_Group,
    Relation_to_Island
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

# View first 5 CpG sites
cat("Details of first 5 SERPINB1-related CpG sites:\n")
print(head(cpgs_serpinb1_info, 5))

# Save CpG site information
write.csv(cpgs_serpinb1_info, "serpinb1_related_cpg_sites.csv", row.names = FALSE)
cat("Saved SERPINB1-related CpG site information to: serpinb1_related_cpg_sites.csv\n")

# --------------------- 5. Integrate Methylation Data with Clinical Data ---------------------
cat("Step 4: Integrating methylation data with clinical data...\n")

# Filter methylation data for SERPINB1-related CpG sites
meth_serpinb1 <- beta_df %>%
  filter(CpG %in% cpg_list)%>%
  select(CpG ,all_of(clinical_df$sample_id)) %>%
  filter(rowSums(!is.na(select(., -CpG))) > 0)
    # Keep only target CpG sites

# Convert to long format
meth_long <- meth_serpinb1 %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample_id",
    values_to = "Beta_value"
  ) %>%
  filter(!is.na(Beta_value))  # Remove missing values

# Save filtered methylation data
write.csv(meth_long, "serpinb1_cpg_methylation_long_format.csv", row.names = FALSE)
cat("Saved filtered methylation data to: serpinb1_cpg_methylation_long_format.csv\n")

# Modified complete merging process
library(dplyr)

# 1. Process methylation data
meth_long <- meth_long %>%
  mutate(
    patient_id = substr(sample_id, 1, 12),  # Extract patient ID
    sample_type = case_when(
      grepl("-01[^-]*$", sample_id) ~ "Primary",
      grepl("-06[^-]*$", sample_id) ~ "Normal",
      grepl("-11[^-]*$", sample_id) ~ "Normal",
      TRUE ~ "Other"
    )
  )

# 2. Process clinical data
clinical_df <- clinical_df %>%
  mutate(patient_id = substr(sample_id, 1, 12))

# 3. Merge data
merged_data <- meth_long %>%
  inner_join(
    clinical_df %>% select(patient_id, Gleason, sample_type_clin = sample_type),
    by = "patient_id"
  ) %>%
  left_join(
    cpgs_serpinb1_info,
    by = c("CpG" = "Name")
  )

# 4. Validate and clean results
if (nrow(merged_data) == 0) {
  # Detailed diagnostics
  cat("Methylation patient IDs:", paste(unique(meth_long$patient_id), collapse = ", "), "\n")
  cat("Clinical patient IDs:", paste(unique(clinical_df$patient_id), collapse = ", "), "\n")
  cat("Intersection:", paste(intersect(unique(meth_long$patient_id), unique(clinical_df$patient_id)), collapse = ", "), "\n")
  
  stop("Merge failed: No common patient IDs")
} else {
  cat("Successfully merged", nrow(merged_data), "rows of data, containing",
      n_distinct(merged_data$patient_id), "patients and",
      n_distinct(merged_data$CpG), "CpG sites\n")
  
  # View results
  print(head(merged_data))
}

# Data validation
cat("Merged data dimensions:", dim(merged_data), "\n")
cat("Number of samples:", length(unique(merged_data$sample_id)), "\n")
cat("Number of CpG sites:", length(unique(merged_data$CpG)), "\n")

# Save merged data
write.csv(merged_data, "merged_methylation_clinical_data.csv", row.names = FALSE)
cat("Saved merged data to: merged_methylation_clinical_data.csv\n")

# --------------------- 6. Calculate Spearman Correlations ---------------------
cat("Step 5: Calculating Spearman correlation between each CpG site and Gleason score...\n")

# 1. Check missing values (key variables)
cat("Beta_value missing values:", sum(is.na(merged_data$Beta_value)), "\n")
cat("Gleason missing values:", sum(is.na(merged_data$Gleason)), "\n")

# 2. Count valid observations per CpG group (both Beta_value and Gleason non-missing)
group_obs <- merged_data %>%
  filter(!is.na(Beta_value), !is.na(Gleason)) %>%  # Keep only valid pairs
  group_by(CpG, Position, Gene_Region) %>%
  summarise(valid_n = n(), .groups = "drop")

# 3. Check groups with insufficient observations (<3, cannot calculate correlation)
low_obs_groups <- group_obs %>% filter(valid_n < 3)
if (nrow(low_obs_groups) > 0) {
  cat("⚠️ Found", nrow(low_obs_groups), "CpG groups with insufficient observations (<3) for correlation calculation:\n")
  print(head(low_obs_groups))
}

cor_results <- merged_data %>%
  group_by(CpG, Position, Gene_Region) %>%  # Group by CpG
  summarise(
    Spearman_corr = cor.test(Beta_value, Gleason, method = "spearman")$estimate,
    P_value = cor.test(Beta_value, Gleason, method = "spearman")$p.value,
    N = n(),  # Sample count
    .groups = "drop"
  ) %>%
  mutate(
    # Add significance markers
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Add correlation direction markers
    Direction = ifelse(Spearman_corr > 0, "Positive", "Negative")
  ) %>%
  arrange(desc(abs(Spearman_corr)))  # Sort by absolute correlation value

# View results
cat("Correlation analysis results summary:\n")
print(head(cor_results))

# Save results to CSV file
write.csv(cor_results, "SERPINB1_Gleason_correlation_results.csv", row.names = FALSE)
cat("Correlation analysis results saved to 'SERPINB1_Gleason_correlation_results.csv'\n")

# --------------------- 7. Visualize Results (Heatmap) ---------------------
cat("Step 6: Generating heatmap...\n")

# Prepare heatmap data
heatmap_data <- cor_results %>%
  select(CpG, Spearman_corr) %>%
  arrange(Spearman_corr)  # Sort for ordered heatmap display

# Convert to matrix format
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$CpG

# Add annotation information
row_anno <- cor_results %>%
  select(CpG, Gene_Region, Significance) %>%
  distinct() %>%
  column_to_rownames("CpG")

# Create gene region annotation vector
gene_region_vec <- row_anno$Gene_Region 
names(gene_region_vec) <- rownames(row_anno)

# Create significance annotation vector
significance_vec <- row_anno$Significance
names(significance_vec) <- rownames(row_anno)

# 1. Handle empty strings in Significance
significance_vec[significance_vec == ""] <- "NS"

# 2. Define color mappings
gene_region_colors <- c(
  "Promoter" = "#E41A1C", 
  "Gene Body" = "#377EB8", 
  "5'UTR" = "#4DAF4A", 
  "3'UTR" = "#984EA3",
  "Other" = "#FF7F00"
)

significance_colors <- c(
  "***" = "black",
  "**" = "darkgray",
  "*" = "lightgray",
  "NS" = "white"  # Not significant
)

# 3. Create row annotations
ha <- rowAnnotation(
  Gene_Region = gene_region_vec,
  Significance = significance_vec,
  col = list(
    Gene_Region = gene_region_colors,
    Significance = significance_colors
  ),
  annotation_name_side = "top",
  annotation_legend_param = list(
    Gene_Region = list(title = "Genomic Region"),
    Significance = list(title = "Significance")
  )
)

# Create heatmap
ht <- Heatmap(
  heatmap_matrix,
  name = "Spearman\nCorrelation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = FALSE,  # Don't cluster rows, maintain sorting
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "Correlation between SERPINB1 CpG Methylation and Gleason Score",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  right_annotation = ha,  # Add right annotation
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal"
  )
)

# Plot and save results
png("SERPINB1_Gleason_correlation_heatmap2.png",  width = 1200,   
    height = 800 , res=100  )
#draw(ht, heatmap_legend_side = "left", annotation_legend_side = "right")
draw(ht, 
     heatmap_legend_side = "left", 
     annotation_legend_side = "right",
     padding = unit(c(2, 10,2, 2), "mm"),
     )
dev.off()

cat("Heatmap saved to 'SERPINB1_Gleason_correlation_heatmap.png'\n")
cat("Analysis pipeline completed successfully!\n")

