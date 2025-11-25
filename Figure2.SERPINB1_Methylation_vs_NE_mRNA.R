# ======================================================================
# TCGA-PRAD Methylation and Expression Correlation Analysis
# Objective: Analyze association between SERPINB1 promoter methylation 
#            and expression of MMP2, MMP9, SNAIL, TWIST, ERK1, ERK2
# Process improvements:
#   1. Fixed data type conversion issues (character → numeric)
#   2. Unified sample ID handling logic
#   3. Enhanced error handling and logging
#   4. Optimized visualization parameters
# ======================================================================

# -------------------------- 1. Load Dependencies --------------------------
suppressPackageStartupMessages({
  library(data.table)    # Efficient reading of large data files (supports .gz compression)
  library(dplyr)         # Data cleaning and transformation
  library(tidyr)         # Data reshaping (long/wide format)
  library(ggplot2)      # Advanced plotting system
  library(cowplot)      # Plot composition and layout
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Illumina 450K methylation array annotation
})
message("All dependencies loaded successfully")
setwd("Figure2_Data")
# -------------------------- 2. Read Methylation Data --------------------------
# 2.1 Define file paths
meth_file <- "data/TCGA.PRAD.sampleMap_HumanMethylation450.gz"

# 2.2 Read methylation data (using data.table for speed)
message("Reading methylation data...")
methylation_data <- data.table::fread(
  cmd = paste("zcat", meth_file),  # Use zcat on Linux/macOS; 'gzip -d -c' on Windows
  header = FALSE,
  sep = "\t",
  data.table = FALSE,
  na.strings = c("NA", "N/A", "", ".")  # Specify possible missing value representations
)

# 2.3 Extract sample IDs and CpG sites with data type conversion
message("Processing methylation data structure...")
# Sample IDs (row 1, columns 2+)
sample_names_meth <- as.character(methylation_data[1, -1])
# CpG sites (column 1, rows 2+)
cpg_names <- as.character(methylation_data[-1, 1])

# 2.4 Build Beta value matrix (ensure numeric type)
beta_matrix <- methylation_data[-1, -1] %>%
  as.matrix() %>%
  apply(2, function(x) {
    # Replace any non-numeric characters with NA
    x <- gsub("[^0-9.]", NA, x)
    as.numeric(x)
  })

# Set row and column names
rownames(beta_matrix) <- cpg_names
colnames(beta_matrix) <- sample_names_meth

# 2.5 Convert to dataframe (add CpG column)
beta_df <- as.data.frame(beta_matrix) %>%
  tibble::rownames_to_column(var = "CpG")

message(sprintf("Successfully read methylation data for %d CpG sites across %d samples", 
                nrow(beta_df), ncol(beta_matrix)))

# -------------------------- 3. Extract SERPINB1 Promoter Regions --------------------------
# 3.1 Get 450K methylation array annotation
message("Loading methylation array annotation data...")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 3.2 Filter for SERPINB1 promoter CpG sites
cpgs_promoter_serpinb1 <- as.data.frame(anno) %>%
  filter(
    grepl("^SERPINB1$|;SERPINB1|^SERPINB1;", UCSC_RefGene_Name),  # Exact match for SERPINB1 gene
    grepl("TSS1500|TSS200|5'UTR|1stExon", UCSC_RefGene_Group)  # Promoter-related regions
  ) %>%
  # Add new column: explicitly label promoter region type
  mutate(
    Promoter_Type = case_when(
      grepl("TSS1500", UCSC_RefGene_Group) ~ "TSS1500",
      grepl("TSS200", UCSC_RefGene_Group) ~ "TSS200",
      grepl("5'UTR", UCSC_RefGene_Group) ~ "5'UTR",
      grepl("1stExon", UCSC_RefGene_Group) ~ "1stExon",
      TRUE ~ "Other"
    )
  )

# 3.3 Get CpG site list
cpg_list_promoter <- unique(cpgs_promoter_serpinb1$Name)

# 3.4 Validate target CpG sites found
if (length(cpg_list_promoter) == 0) {
  stop("Error: No CpG sites found in SERPINB1 promoter regions! Please check annotation data.")
} else {
  message(sprintf("Found %d SERPINB1 promoter-associated CpG sites", length(cpg_list_promoter)))
  message("Promoter region type distribution:")
  print(table(cpgs_promoter_serpinb1$Promoter_Type))
}

# 3.5 Extract target CpG Beta values and calculate average methylation
message("Calculating average SERPINB1 promoter methylation level...")
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
    CpG_count = n(),
    Non_NA_count = sum(!is.na(Beta_value))
  ) %>%
  # Filter out samples where all CpG sites are NA
  filter(Non_NA_count > 0)

message(sprintf("Calculated average SERPINB1 promoter methylation for %d samples (filtered out %d invalid samples)",
                nrow(avg_meth_promoter), length(cpg_list_promoter) - nrow(avg_meth_promoter)))

# -------------------------- 4. Read mRNA Expression Data --------------------------
# 4.1 Define file path
expr_file <- "data/TCGA.PRAD.sampleMap_HiSeqV2.gz"

# 4.2 Read expression data
message("Reading expression data...")
expr_data <- data.table::fread(
  cmd = paste("zcat", expr_file),
  header = FALSE,
  sep = "\t",
  data.table = FALSE,
  na.strings = c("NA", "N/A", "", "null")
)

# 4.3 Process gene and sample names
message("Processing expression data structure...")
gene_names <- as.character(expr_data[-1, 1])  # Remove header row
sample_names <- as.character(expr_data[1, -1])

# 4.4 Build expression matrix (convert to numeric)
expr_matrix <- expr_data[-1, -1] %>%
  as.matrix() %>%
  apply(2, function(x) {
    suppressWarnings(as.numeric(x))  # Suppress potential warnings
  })

# Set row and column names
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- sample_names

# 4.5 Extract target gene expression
target_genes <- c("SERPINB1", "MMP2", "MMP9",  "SNAI1", "TWIST1","MAPK3", "MAPK1") # MAPK3=ERK1, MAPK1=ERK2
expr_target <- as.data.frame(expr_matrix) %>%
  tibble::rownames_to_column(var = "Gene") %>%
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

message(sprintf("Successfully extracted expression data for %d samples (genes: %s)", 
                nrow(expr_target), paste(target_genes, collapse=", ")))

# -------------------------- 5. Read Clinical Data --------------------------
# 5.1 Define file path and key column names
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"
sample_id_col <- "sampleID"      # Sample ID column name
sample_type_col <- "sample_type" # Sample type column name

# 5.2 Read clinical data
message("Reading clinical data...")
clinical_data <- data.table::fread(
  file = clinical_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE,
  fill = TRUE  # Handle potentially irregular rows
)

# 5.3 Preprocess clinical data
message("Processing clinical data...")
clinical_df <- clinical_data %>%
  # Select key columns and rename
  select(
    sample_id = !!sym(sample_id_col),
    sample_type = !!sym(sample_type_col)
  ) %>%
  # Remove duplicate samples
  distinct(sample_id, .keep_all = TRUE) %>%
  # Create grouping variable
  mutate(
    Group = case_when(
      sample_type == "Solid Tissue Normal" ~ "Normal tissue",
      grepl("Tumor|Primary", sample_type, ignore.case = TRUE) ~ "Cancer tissue",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter invalid groups and missing values
  filter(!is.na(Group))

# 5.4 Check group distribution
group_counts <- table(clinical_df$Group)
message("\nClinical data group distribution:")
message(sprintf(" - Cancer tissue: %d samples", group_counts["Cancer tissue"]))
message(sprintf(" - Normal tissue: %d samples", group_counts["Normal tissue"]))

# -------------------------- 6. Merge All Three Datasets --------------------------
# 6.1 Ensure consistent sample ID format
message("\nMerging datasets...")
# Extract all sample IDs
all_samples <- list(
  methylation = unique(avg_meth_promoter$sample_id),
  expression = unique(expr_target$sample_id),
  clinical = unique(clinical_df$sample_id)
)

# 6.2 Find common sample IDs
common_samples <- Reduce(intersect, all_samples)

# 6.3 Validate common samples exist
if (length(common_samples) == 0) {
  # Detailed error diagnostics
  warning("Methylation sample examples: ", paste(head(all_samples$methylation, 3), collapse = ", "))
  warning("Expression sample examples: ", paste(head(all_samples$expression, 3), collapse = ", "))
  warning("Clinical sample examples: ", paste(head(all_samples$clinical, 3), collapse = ", "))
  stop("Error: No common sample IDs found across methylation, expression and clinical data! Check data formats.")
} else {
  m<- sprintf("\nFound %d common samples for downstream analysis", length(common_samples))
  message(m)
}

# 6.4 Merge data
merged_data <- avg_meth_promoter %>%
  inner_join(expr_target, by = "sample_id") %>%
  inner_join(clinical_df, by = "sample_id")

# Check final dataset
na_stats <- colSums(is.na(merged_data))
message(sprintf("Final dataset contains %d samples", nrow(merged_data)))
message("\nMissing value statistics:")
print(na_stats)


# -------------------------- 7. Correlation Analysis and Visualization --------------------------

# -------------------------- 7. Correlation Analysis (Overall) --------------------------

# 7.1 Enhanced visualization function (non-grouped correlation)
plot_overall_cor_scatter <- function(
    data, 
    x_var, 
    y_var,
    color_var = "Group",  # Still color by group but calculate overall correlation
    color_palette = c("Normal tissue" = "#377EB8", "Cancer tissue" = "#D55E00"),
    title = NULL,
    x_lab = "SERPINB1 promoter methylation (β-value)",
    y_lab = NULL
) {
  
  # Validate required columns exist
  required_cols <- c(x_var, y_var, color_var)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create filtered plot data
  plot_data <- data %>%
    select(all_of(c(x_var, y_var, color_var))) %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  
  if (nrow(plot_data) == 0) {
    warning(sprintf("No complete data available for plotting %s vs %s", y_var, x_var))
    return(NULL)
  }
  
  # Calculate overall correlation statistics (not by group)
  corr_result <- plot_data %>%
    summarise(
      Pearson_r = cor(!!sym(x_var), !!sym(y_var), method = "pearson", use = "complete.obs"),
      P_value = suppressWarnings(cor.test(!!sym(x_var), !!sym(y_var), method = "pearson")$p.value),
      Sample_size = n()
    ) %>%
    mutate(
      Pearson_r = round(Pearson_r, 3),
      P_value = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)),
      # Calculate annotation position (top-right corner)
      x_pos = max(plot_data[[x_var]], na.rm = TRUE) * 0.9,
      y_pos = max(plot_data[[y_var]], na.rm = TRUE) * 0.9
    )
  
  # Handle titles and labels
  plot_title <- ifelse(is.null(title), 
                       paste(y_var, "expression vs SERPINB1 methylation"),
                       title)
  
  y_lab <- ifelse(is.null(y_lab), 
                  paste(y_var, "mRNA expression"),  #paste(y_var, "mRNA expression (log2(TPM+1))"),
                  y_lab)
  
  # Create scatter plot (still colored by group)
  p <- ggplot(plot_data, aes(
    x = .data[[x_var]], 
    y = .data[[y_var]], 
    color = .data[[color_var]]
  )) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(
      method = "lm", 
      se = TRUE, 
      color = "black",  # Single trend line for overall correlation
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_color_manual(values = color_palette) +
    labs(
      title = plot_title,
      x = x_lab,
      y = y_lab,
      color = "Tissue type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  # Add single statistical annotation (overall correlation)
  p <- p + annotate(
    "text",
    x = corr_result$x_pos,
    y = corr_result$y_pos,
    label = sprintf("Overall:\nr = %s\nP = %s\nN = %d", 
                    corr_result$Pearson_r, 
                    corr_result$P_value,
                    corr_result$Sample_size),
    color = "black",  # Black for overall correlation
    hjust = 1,
    vjust = 1,
    size = 4,
    fontface = "bold"
  )
  
  # Output statistical results
  message(sprintf("\n%s vs %s overall correlation statistics:", y_var, x_var))
  print(corr_result)
  
  return(list(plot = p, corr_stats = corr_result))
}

# 7.2 Generate non-grouped correlation plots for all target genes
message("\nGenerating overall correlation plots...")

# Create list to store plots
cor_plots_overall <- list()

# Generate plots for each target gene
for (gene in target_genes[-1]) { # Skip SERPINB1 (x-axis variable)
  p <- plot_overall_cor_scatter(
    data = merged_data, 
    x_var = "SERPINB1_Promoter_Beta", 
    y_var = gene,
    title = paste(gene, "expression vs SERPINB1 methylation")
  )
  
  if (!is.null(p)) {
    cor_plots_overall[[gene]] <- p
  }
}

# 7.3 Combine and save non-grouped correlation plots
if (length(cor_plots_overall) > 0) {
  # Create plot grid with consistent sizing
  plot_grid <- cowplot::plot_grid(
    plotlist = lapply(cor_plots_overall, function(x) x$plot),
    ncol = 2,
    labels = paste0(LETTERS[1:length(cor_plots_overall)], ": ", names(cor_plots_overall)),
    label_size = 12,
    align = "hv",
    axis = "lb"  # Align axes for clean layout
  )
  
  # Add title and save with proper dimensions
  title <- ggdraw() + 
    draw_label("SERPINB1 Promoter Methylation Correlations (Overall)", 
               size = 14, fontface = "bold", x = 0.5, hjust = 0.5)
  
  final_plot <- plot_grid(
    title,
    plot_grid,
    ncol = 1,
    rel_heights = c(0.1, 1)  # Title takes 10% height
  )
  
  # Calculate dynamic height based on number of plots
  n_plots <- length(cor_plots_overall)
  plot_rows <- ceiling(n_plots/2)
  plot_height <- max(5, 3 * plot_rows)  # Minimum 5 inches, 3 inches per row
  
  # Save high-quality PDF
  ggsave(
    filename = "SERPINB1_Methylation_Correlations_Overall.pdf",
    plot = final_plot,
    width = 14,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
  
  # Save PNG for quick viewing
  ggsave(
    filename = "SERPINB1_Methylation_Correlations_Overall.png",
    plot = final_plot,
    width = 14,
    height = plot_height,
    dpi = 150
  )
  
  message(sprintf(
    "\nSuccess! Saved combined plot with %d subplots to:\n- SERPINB1_Methylation_Correlations_Overall.pdf\n- SERPINB1_Methylation_Correlations_Overall.png",
    n_plots
  ))
  
  # Save overall correlation results
  overall_results <- lapply(cor_plots_overall, function(x) {
    x$corr_stats %>% mutate(Gene = x$corr_stats$y_var[1])
  }) %>% bind_rows()
  
  write.csv(overall_results, "SERPINB1_overall_correlation_results.csv", row.names = FALSE)
  message("Saved overall correlation results to: SERPINB1_overall_correlation_results.csv")
  
} else {
  warning("No valid plots generated in cor_plots_overall list")
}


######################
# -------------------------- 2. Save Methylation Data --------------------------
# After beta_df creation (2.5)
write.csv(beta_df, "methylation_beta_values_all_cpgs.csv", row.names = FALSE)
message("Saved all CpG beta values to: methylation_beta_values_all_cpgs.csv")

# Save SERPINB1 promoter CpGs info
write.csv(cpgs_promoter_serpinb1, "serpinb1_promoter_cpgs_annotation.csv", row.names = FALSE)
message("Saved SERPINB1 promoter CpG annotations to: serpinb1_promoter_cpgs_annotation.csv")

# Save average promoter methylation
write.csv(avg_meth_promoter, "serpinb1_promoter_avg_methylation.csv", row.names = FALSE)
message("Saved SERPINB1 average promoter methylation to: serpinb1_promoter_avg_methylation.csv")

# -------------------------- 4. Save Expression Data --------------------------
# After expr_matrix creation (4.4)
write.csv(expr_matrix, "expression_matrix_all_genes.csv", row.names = TRUE)
message("Saved full expression matrix to: expression_matrix_all_genes.csv")

# After expr_target creation (4.5)
write.csv(expr_target, "target_genes_expression_data.csv", row.names = FALSE)
message("Saved target genes expression data to: target_genes_expression_data.csv")

# -------------------------- 5. Save Clinical Data --------------------------
# After clinical_df creation (5.3)
write.csv(clinical_df, "clinical_data_with_groups.csv", row.names = FALSE)
message("Saved clinical data with groups to: clinical_data_with_groups.csv")

# -------------------------- 6. Save Merged Data --------------------------
# Before merging (sample lists)
sample_lists <- list(
  methylation = data.frame(sample_id = all_samples$methylation),
  expression = data.frame(sample_id = all_samples$expression),
  clinical = data.frame(sample_id = all_samples$clinical)
)
write.csv(sample_lists$methylation, "methylation_sample_list.csv", row.names = FALSE)
write.csv(sample_lists$expression, "expression_sample_list.csv", row.names = FALSE)
write.csv(sample_lists$clinical, "clinical_sample_list.csv", row.names = FALSE)
message("Saved sample ID lists for each data type")

# Save common samples
write.csv(data.frame(sample_id = common_samples), "common_samples_across_datasets.csv", row.names = FALSE)
message("Saved common sample IDs to: common_samples_across_datasets.csv")

# After merging (6.4)
write.csv(merged_data, "merged_methylation_expression_clinical_data.csv", row.names = FALSE)
message("Saved final merged dataset to: merged_methylation_expression_clinical_data.csv")

# -------------------------- 7. Save Correlation Results --------------------------
# After correlation analysis (7.2)
correlation_results <- lapply(cor_plots, function(p) {
  if (!is.null(p$corr_stats)) {
    return(p$corr_stats)
  }
}) %>% bind_rows(.id = "Gene")

write.csv(correlation_results, "gene_expression_correlation_results.csv", row.names = FALSE)
message("Saved correlation results to: gene_expression_correlation_results.csv")

# Save plot data for each gene
for (gene in names(cor_plots)) {
  plot_data <- merged_data %>%
    select(sample_id, Group, SERPINB1_Promoter_Beta, !!sym(gene)) %>%
    rename(Expression = !!sym(gene))
  
  write.csv(plot_data, paste0("plot_data_", gene, "_vs_SERPINB1.csv"), row.names = FALSE)
}
message("Saved individual plot data for each gene comparison")

# -------------------------- 8. Save Analysis Metadata --------------------------
analysis_metadata <- data.frame(
  AnalysisDate = Sys.Date(),
  TCGA_DataVersion = "PRAD HiSeqV2",
  MethylationPlatform = "Illumina HumanMethylation450",
  TargetGenes = paste(target_genes, collapse = ", "),
  NormalSamples = sum(merged_data$Group == "Normal tissue"),
  CancerSamples = sum(merged_data$Group == "Cancer tissue"),
  SERPINB1_CpGs = length(cpg_list_promoter)
)

write.csv(analysis_metadata, "analysis_metadata.csv", row.names = FALSE)
message("Saved analysis metadata to: analysis_metadata.csv")

# -------------------------- 9. Create Consolidated Output --------------------------
# Create timestamped directory for outputs
output_dir <- paste0("TCGA_PRAD_SERPINB1_Analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir)

# Move all CSV files to output directory
csv_files <- list.files(pattern = "\\.csv$")
file.copy(csv_files, output_dir)
file.remove(csv_files)  # Clean up working directory

message(paste("\nAll intermediate files saved to directory:", output_dir))
message("Files created:")
message(paste("-", list.files(output_dir), collapse = "\n"))

# [Rest of the original plotting code remains unchanged]

