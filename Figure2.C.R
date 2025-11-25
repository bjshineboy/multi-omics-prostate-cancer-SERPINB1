# Load required packages (with purpose descriptions)
library(data.table)    # Efficient reading of large files (supports .gz format)
library(dplyr)         # Data cleaning and transformation
library(tidyr)         # Long-to-wide and wide-to-long format conversion
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Provides CpG site annotations for 450k array (hg19 genome)
setwd("Figure2_Data")
# Define methylation data path (clarify data source)
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

# Read methylation data (.gz format requires zcat decompression, cmd parameter enables pipe reading)
methylation_data <- data.table::fread(
  cmd = paste("zcat", meth_file),  # Use system command zcat to decompress .gz files
  header = TRUE,                   # First row contains column names (sample IDs)
  sep = "\t",                      # Tab delimiter
  data.table = FALSE               # Return data.frame format (better for dplyr operations)
)

# Save raw methylation data
write.csv(methylation_data, "raw_methylation_data.csv", row.names = FALSE)
cat("Saved raw methylation data to: raw_methylation_data.csv\n")

# Extract sample IDs (methylation data columns 2+ in row 1 are sample IDs, column 1 is "CpG" identifier)
sample_names_meth <- colnames(methylation_data)[-1]
# Extract CpG site names (column 1 rows 2+ are CpG IDs, row 1 is header)
cpg_names <- as.character(methylation_data[-1, 1])

# Build Beta value matrix (core methylation data): rows=CpG sites, columns=samples, values=methylation Beta values (0=unmethylated, 1=fully methylated)
beta_matrix <- as.matrix(methylation_data[-1, -1])  # Exclude row 1 (sample names) and column 1 (CpG names)
rownames(beta_matrix) <- cpg_names                  # Row names as CpG IDs
colnames(beta_matrix) <- sample_names_meth          # Column names as sample IDs

# ---------------------  Read Clinical Data  ---------------------
# Clinical data file path
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"

# Read clinical data
clinical_data <-  fread(
      clinical_file,
      header = TRUE,
      sep = "\t",
      data.table = FALSE,
      na.strings = c("", "NA", "[Not Available]")  # Specify missing value representations
    )
 
clinical_df <- clinical_data %>%
  select(sampleID, sample_type) %>%
  filter(sample_type=="Primary Tumor")

# ------------------------------------------ end  -----------

# Convert to data.frame (for easier filtering), keeping CpG ID as a column
beta_df <- as.data.frame(beta_matrix) %>%
  select(all_of(clinical_df$sampleID)) %>% 
  filter(rowSums(!is.na(select(., -CpG))) > 0) %>%
  tibble::rownames_to_column(var = "CpG")  # Add "CpG" column (original row names)

# Save beta matrix
#write.csv(beta_df, "methylation_beta_matrix.csv", row.names = FALSE)
cat("Saved methylation beta matrix to: methylation_beta_matrix.csv\n")

# Extract CpG sites associated with SERPINB1 gene from 450k annotation library (key step: link genes to CpGs)
# Annotation library contains genomic positions, annotated genes etc. for each CpG site
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Get annotation data
cpgs_serpinb1 <- as.data.frame(anno) %>%
  filter(grepl("SERPINB1\\b", UCSC_RefGene_Name))  # Filter CpG sites annotated with "SERPINB1" (UCSC_RefGene_Name contains gene annotations)

# Save SERPINB1 CpG annotations
write.csv(cpgs_serpinb1, "serpinb1_cpg_annotations.csv", row.names = FALSE)
cat("Saved SERPINB1 CpG annotations to: serpinb1_cpg_annotations.csv\n")

# Extract list of CpG IDs for SERPINB1 gene (for filtering methylation data)
cpg_list <- cpgs_serpinb1$Name  # "Name" column contains CpG IDs (e.g. "cg00012345")
# Filter methylation data to keep only SERPINB1 gene CpG sites
meth_serpinb1 <- beta_df %>%
  filter(CpG %in% cpg_list)  # Keep rows where "CpG" column is in cpg_list

# Save filtered methylation data
write.csv(meth_serpinb1, "serpinb1_cpg_methylation.csv", row.names = FALSE)
cat("Saved SERPINB1 CpG methylation data to: serpinb1_cpg_methylation.csv\n")


# Convert from wide (rows=CpGs, columns=samples) to long format (rows=sample-CpG pairs) for calculating average methylation per sample
meth_long <- meth_serpinb1 %>%
  pivot_longer(
    cols = -CpG,               # Keep "CpG" column, convert other columns (samples) to rows
    names_to = "sample_id",    # New column "sample_id" stores original column names (sample IDs)
    values_to = "Beta_value"   # New column "Beta_value" stores original cell values (methylation Beta values)
  )

# Calculate average methylation level per sample for SERPINB1 gene (mean of multiple CpG sites represents overall gene methylation state)
avg_meth <- meth_long %>%
  group_by(sample_id) %>%          # Group by sample
  summarise(
    SERPINB1_Beta = mean(Beta_value, na.rm = TRUE)  # Calculate mean Beta value for all SERPINB1 CpGs per sample (ignoring NAs)
  )

# Save average methylation
write.csv(avg_meth, "serpinb1_average_methylation.csv", row.names = FALSE)
cat("Saved SERPINB1 average methylation to: serpinb1_average_methylation.csv\n")

# Output key statistics (validate data processing)
cat("Successfully extracted", length(cpg_list), "CpG sites for SERPINB1 gene, calculated average methylation Beta values for", nrow(avg_meth), "samples\n")

# Define expression data path (TCGA HiSeqV2.gz file containing gene expression)
expr_file <- "data/TCGA.PRAD.sampleMap_HiSeqV2.gz"

# Read expression data (similar processing logic as methylation data)
expr_data <- data.table::fread(
  cmd = paste("zcat", expr_file),  # Use zcat to decompress .gz file
  header = TRUE,                   # First row contains column names (sample IDs)
  sep = "\t",                      # Tab delimiter
  data.table = FALSE               # Return data.frame format
)

# Save raw expression data
write.csv(expr_data, "raw_expression_data.csv", row.names = FALSE)
cat("Saved raw expression data to: raw_expression_data.csv\n")

# Extract gene names (expression data column 1 contains gene symbols like "TP53")
gene_names_expr <- as.character(expr_data[, 1])  # Column 1 is gene names
# Extract sample IDs (expression data row 1 columns 2+ are sample IDs)
sample_names_expr <- colnames(expr_data)[-1]

# Build expression matrix (rows=genes, columns=samples, values=expression levels (log-transformed FPKM/TPM))
expr_matrix <- as.matrix(expr_data[, -1])  # Exclude row 1 (sample names) and column 1 (gene names)
rownames(expr_matrix) <- gene_names_expr     # Row names as gene symbols
colnames(expr_matrix) <- sample_names_expr   # Column names as sample IDs

# Save expression matrix
write.csv(expr_matrix, "expression_matrix.csv", row.names = TRUE)
cat("Saved expression matrix to: expression_matrix.csv\n")

# Convert to data.frame (for easier filtering)
expr_df <- as.data.frame(expr_matrix)

# Define target genes for NE-MAPK signaling axis (based on research hypothesis, these genes are related to prostate cancer progression)
# Note: MMP2/9=matrix metalloproteinases (extracellular matrix degradation); SNAIL/TWIST=EMT transcription factors (promote metastasis); MAPK3/1=ERK1/2 (core kinases in MAPK pathway)
target_genes <- c("MMP2", "MMP9", "SNAI1", "TWIST1", "MAPK3", "MAPK1")

# Filter expression data for target genes (keep genes of interest)
expr_target <- expr_df %>%
  tibble::rownames_to_column(var = "Gene") %>%  # Add "Gene" column (original row names=gene symbols)
  filter(Gene %in% target_genes) %>%     # Keep rows where Gene is in target_genes
  pivot_longer(
    cols = -Gene,               # Keep "Gene" column, convert other columns (samples) to rows
    names_to = "sample_id",    # New column "sample_id" stores original column names (sample IDs)
    values_to = "Expression"   # New column "Expression" stores original cell values (expression levels)
  )

# Save target gene expression
write.csv(expr_target, "target_gene_expression_long.csv", row.names = FALSE)
cat("Saved target gene expression (long format) to: target_gene_expression_long.csv\n")

# Convert from long to wide format (one row per sample, one column per gene) for merging with methylation data
expr_wide <- expr_target %>%
  pivot_wider(
    names_from = Gene,           # Use "Gene" column values as new column names (e.g. "MMP2")
    values_from = Expression    # Use "Expression" column values to fill new columns
  )

# Save wide format expression
write.csv(expr_wide, "target_gene_expression_wide.csv", row.names = FALSE)
cat("Saved target gene expression (wide format) to: target_gene_expression_wide.csv\n")

# Key check: match sample IDs between methylation and expression data (correlation analysis requires consistent samples)
common_samples <- intersect(avg_meth$sample_id, expr_wide$sample_id)  # Get intersection of sample IDs  **** Here miss one sample

# If no common samples, terminate program (indicates sample ID format mismatch, need to check data sources)
if (length(common_samples) == 0) {
  stop("Methylation and expression data samples don't match at all! Please check if sample IDs follow TCGA standard format (e.g. 'TCGA-XX-XXXX-01A').")
}

# Filter methylation data to keep only common samples
avg_meth_common <- avg_meth %>%
  filter(sample_id %in% common_samples)  # Keep rows where sample_id is in common_samples

# Filter expression data to keep only common samples
expr_common <- expr_wide %>%
  filter(sample_id %in% common_samples)

# Save filtered data
write.csv(avg_meth_common, "serpinb1_methylation_common_samples.csv", row.names = FALSE)
write.csv(expr_common, "target_gene_expression_common_samples.csv", row.names = FALSE)
cat("Saved filtered data for common samples\n")

# Output number of common samples (validate matching results)
cat("Number of common samples between methylation and expression data:", nrow(avg_meth_common), "(needs to be large enough for reliable correlations)\n")

# Merge methylation data (avg_meth_common) and expression data (expr_common): join by sample_id
# Result: one row per sample, containing SERPINB1_Beta (methylation) and target gene expression
merged_data <- inner_join(
  avg_meth_common,  # Methylation data (sample_id + SERPINB1_Beta)
  expr_common,      # Expression data (sample_id + target gene expression)
  by = "sample_id"  # Join key: sample_id (must match)
)

# Save merged data
write.csv(merged_data, "merged_methylation_expression_data.csv", row.names = FALSE)
cat("Saved merged methylation and expression data to: merged_methylation_expression_data.csv\n")

# Convert from wide to long format (one row per sample-gene pair) for calculating correlations by gene
long_data <- merged_data %>%
  pivot_longer(
    cols = -c(sample_id, SERPINB1_Beta),  # Keep sample_id and SERPINB1_Beta, convert other columns (target genes) to rows
    names_to = "Gene",                    # New column "Gene" stores original column names (e.g. "MMP2")
    values_to = "Expression"              # New column "Expression" stores original cell values (expression)
  )

# Save long format data
write.csv(long_data, "merged_data_long_format.csv", row.names = FALSE)
cat("Saved merged data in long format to: merged_data_long_format.csv\n")

# Calculate Spearman correlation between each gene and SERPINB1 methylation (non-parametric, suitable for non-normal distributions)
cor_results <- long_data %>%
  group_by(Gene) %>%  # Group by gene (calculate correlation separately for each gene)
  summarise(
    Spearman_corr = cor(SERPINB1_Beta, Expression, method = "spearman", use = "complete.obs"),  # Spearman correlation (-1 to 1, larger absolute value = stronger correlation)
    P_value = cor.test(SERPINB1_Beta, Expression, method = "spearman", use = "complete.obs")$p.value,  # Significance p-value (<0.05 = statistically significant)
    Sample_size = n()  # Sample size used in calculation (after excluding NAs)
  ) %>%
  arrange(desc(abs(Spearman_corr)))  # Sort by absolute correlation value (descending, to easily see strongest correlations)

# Save correlation results
write.csv(cor_results, "serpinb1_gene_expression_correlations.csv", row.names = FALSE)
cat("Saved correlation results to: serpinb1_gene_expression_correlations.csv\n")

# Output correlation results (clearly show correlation strength and significance for each gene)
cat("Spearman correlation results between SERPINB1 methylation (SERPINB1_Beta) and NE-MAPK gene expression:\n")
print(cor_results, digits = 3)  # Keep 3 decimal places for readability

# Load visualization packages
library(ComplexHeatmap)  # For complex heatmaps (supports clustering, color customization)
library(circlize)       # For color gradient settings

# Convert correlation results to matrix (rows=genes, columns=correlation metrics, values=Spearman_corr)
# Note: Here we only keep Spearman_corr as heatmap values (can add p-values if needed)
heatmap_mat <- cor_results %>%
  select(Gene, Spearman_corr) %>%
  pivot_wider(names_from = Gene, values_from = Spearman_corr) %>%
  as.matrix()  # Convert to matrix (rows=1 row, columns=genes, values=correlation coefficients)

# Custom color gradient (red=negative correlation, white=no correlation, blue=positive correlation)
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),  # Value range for colors (-1 to 1)
  colors = c("red", "white", "blue")  # Negative=red, positive=blue
)

# Draw heatmap (clearly shows correlation patterns between genes and methylation)
ht_cor <- Heatmap(
  heatmap_mat,               # Input matrix (rows=metrics, columns=genes)
  name = "Spearman\nCorrelation",  # Color legend title (line break improves readability)
  col = col_fun,             # Use custom color gradient
  cluster_rows = FALSE,      # Don't cluster rows (only 1 row: Spearman_corr)
  cluster_columns = TRUE,    # Cluster columns (genes) by correlation strength
  show_row_names = TRUE,     # Show row names ("Spearman_corr")
  show_column_names = TRUE,  # Show column names (gene symbols)
  column_names_gp = gpar(fontsize = 10),  # Column name font size (avoid overlap)
  row_names_gp = gpar(fontsize = 10),     # Row name font size
  cell_fun = function(j, i, x, y, width, height, fill) {  # Custom cell content (add p-values)
    # Extract p-value for current gene (j=column index corresponding to gene)
    p_val <- cor_results$P_value[j]
    # Format p-value (<0.001="***", <0.01="**", <0.05="*")
    p_label <- ifelse(p_val < 0.001, "***", 
                      ifelse(p_val < 0.01, "**", 
                             ifelse(p_val < 0.05, "*", "")))
    # Add p-value label to cell (white text, centered)
    grid.text(p_label, x, y, gp = gpar(col = "white", fontsize = 8))
  },
  column_title = "NE-MAPK Signal Axis Genes",  # Column title (pathway for genes)
  row_title = "SERPINB1 Methylation",         # Row title (methylation metric)
  column_title_side = "top",                  # Column title position (top)
  row_title_side = "left"                     # Row title position (left)
)


png(
  filename = "SERPINB1_methylation_vs_NE-MAPK_expression_correlation_heatmap22.png",
  width = 800,   
  height = 1200,  
  res = 300    
)

 
draw(ht_cor)
 
dev.off()

 
cat("Saved correlation heatmap to: SERPINB1_methylation_vs_NE-MAPK_expression_correlation_heatmap.pdf\n")











