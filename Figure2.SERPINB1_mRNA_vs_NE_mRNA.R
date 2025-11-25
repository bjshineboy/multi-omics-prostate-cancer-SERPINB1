# ====================== Gene Expression Correlation Analysis (TCGA-PRAD) ======================
# Objective: Analyze expression correlations between SERPINB1 and MMP2/MMP9 (stratified by cancer/normal groups)
# Data sources: TCGA-PRAD mRNA expression data (HiSeqV2), clinical sample type data
# Pipeline: Data loading → preprocessing → merging → correlation analysis → visualization
# ==============================================================================


# -------------------------- 1. Load Required Packages --------------------------
# Note: Using conflict-free loading approach
library(data.table, warn.conflicts = FALSE)  # Efficient large file reading
library(dplyr, warn.conflicts = FALSE)       # Data cleaning
library(tidyr, warn.conflicts = FALSE)       # Data reshaping
library(ggplot2)                            # Plotting
library(cowplot)                            # Multi-panel layouts
library(broom)                              # Tidy statistical results

setwd("Figure2_Data")
# -------------------------- 2. Load and Preprocess Expression Data --------------------------
# 2.1 Define file paths and target genes
expr_file <- "data/TCGA.PRAD.sampleMap_HiSeqV2.gz"
target_genes <- c("SERPINB1", "MMP2", "MMP9", "SNAI1", "TWIST1", "MAPK3", "MAPK1")  # MAPK3=ERK1, MAPK1=ERK2

# 2.2 Read data (supports .gz compression)
expr_data <- data.table::fread(
  cmd = paste("zcat", expr_file),  # Use 'gzip -d -c' on Windows
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

# 2.3 Process gene and sample names (TCGA data format adaptation)
gene_names <- expr_data[-1, 1] %>% as.character()  # Column 1: gene names (exclude header)
sample_names <- expr_data[1, -1] %>% as.character()  # Row 1: sample IDs (exclude gene column)

# 2.4 Build expression matrix (rows=genes, columns=samples)
expr_matrix <- as.matrix(expr_data[-1, -1])  # Exclude header and gene column
rownames(expr_matrix) <- gene_names          # Row names: gene names
colnames(expr_matrix) <- sample_names        # Column names: sample IDs

# 2.5 Convert to tidy format (for easier processing)
expr_df <- as.data.frame(expr_matrix) %>%
  tibble::rownames_to_column("Gene") %>%  # Convert gene names to column
  filter(Gene %in% target_genes) %>%      # Keep only target genes
  pivot_longer(                          # Wide to long format (sample-gene-expression)
    cols = -Gene,
    names_to = "sample_id",
    values_to = "Expression"
  ) %>%
  pivot_wider(                           # Long to wide format (sample-gene columns)
    names_from = Gene,
    values_from = Expression
  )

# 2.6 Validate target gene presence
missing_genes <- setdiff(target_genes, colnames(expr_df)[-1])  # Check for missing genes
if (length(missing_genes) > 0) {
  stop(paste("Error: Target genes not found:", paste(missing_genes, collapse = ", ")))
}
message(paste("Successfully extracted expression data for", nrow(expr_df), 
              "samples (genes:", paste(target_genes, collapse = ", "), ")"))


# -------------------------- 3. Load and Preprocess Clinical Data --------------------------
# 3.1 Define file paths and key columns
clinical_file <- "data/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"
sample_id_col <- "sampleID"  # TCGA sample ID column (format: TCGA-XX-XXXX-01)
sample_type_col <- "sample_type"        # Sample type column (e.g., "Solid Tissue Normal")

# 3.2 Read data
clinical_data <- data.table::fread(
  file = clinical_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

# 3.3 Preprocessing: Standardize sample IDs + grouping
clinical_df <- clinical_data %>%
  # Standardize sample IDs (keep first 12 characters to match TCGA format)
  mutate(sample_id = !!sym(sample_id_col)) %>%
  # Select key columns and group (TCGA standard)
  dplyr::select(sample_id, !!sym(sample_type_col)) %>%
  rename(sample_type = !!sym(sample_type_col)) %>%
  mutate(
    Group = case_when(
      sample_type == "Solid Tissue Normal" ~ "Normal",  # Normal tissue (adjacent)
      grepl("Tumor|Primary", sample_type, ignore.case = TRUE) ~ "Cancer",  # Cancer tissue
      TRUE ~ NA_character_  # Other types (e.g., metastases) marked as NA
    )
  ) %>%
  dplyr::filter(!is.na(Group))  # Filter invalid groups

# 3.4 Verify grouping results
message("\nClinical data group distribution:")
print(table(clinical_df$Group))


# -------------------------- 4. Merge Expression and Clinical Data --------------------------
# Merge data and immediately convert target genes to numeric
merged_data <- expr_df %>%
  inner_join(clinical_df, by = "sample_id") %>%
  # Batch convert target genes to numeric (critical step!)
  mutate(across(all_of(target_genes), as.numeric)) %>%
  na.omit()  # Remove samples with NA (optional)

# Verify data types
sapply(merged_data[target_genes], class)

# 4.1 Verify merge results
if (nrow(merged_data) == 0) {
  stop("Error: No data after merging, please check sample ID matching!")
}
message("\nFinal sample count for analysis:", nrow(merged_data))
message("Group distribution after merging:")
print(table(merged_data$Group))

write.csv(merged_data, "")
# -------------------------- 5. Correlation Analysis and Visualization --------------------------
# 5.1 Define plotting function (reusable logic)
plot_cor_scatter <- function(data, x_gene, y_gene, group_var = "Group",
                             x_label = paste(x_gene, "expression (log2(TPM+1))"),
                             y_label = paste(y_gene, "expression (log2(TPM+1))"),
                             title = paste(y_gene, "vs", x_gene, "expression correlation"),
                             color_palette = c("Normal" = "#0072B2", "Cancer" = "#D55E00")) {
  
  # Validate input columns exist (early error catching)
  required_cols <- c(x_gene, y_gene, group_var)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Calculate Pearson correlation (structured results using broom::tidy)
  corr_results <- data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      n = sum(!is.na(.[[x_gene]]) & !is.na(.[[y_gene]])),  # Calculate valid sample size
      cor_test = list(cor.test(.[[x_gene]], .[[y_gene]], method = "pearson"))
    ) %>%
    mutate(
      tidy = map(cor_test, broom::tidy)
    ) %>%
    unnest(tidy) %>%
    dplyr::select(
      Group = !!sym(group_var),
      Pearson_r = estimate,
      P_value = p.value,
      Sample_size = n
    ) %>%
    mutate(
      Pearson_r = round(Pearson_r, 2),
      P_value = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
    )
  
  # Create scatter plot (with confidence interval trend line)
  p <- ggplot(data, aes(
    x = .data[[x_gene]],  # Correct: Use .data[[...]] to avoid variable name conflicts
    y = .data[[y_gene]],
    color = .data[[group_var]]  # Color by group variable (e.g., "Normal"/"Cancer")
  )) +
    geom_point(alpha = 0.6, size = 2) +  # Semi-transparent points (avoid overplotting)
    geom_smooth(
      method = "lm",       # Linear model (basis for Pearson correlation)
      se = TRUE,           # Show 95% confidence interval for trend line
      color = "black",     # Trend line color (black, distinct from points)
      linetype = "dashed"  # Trend line type (dashed for clarity)
    ) +
    # Color scheme (academic convention: normal=blue, cancer=red)
    scale_color_manual(
      values = color_palette,
      name = "Group"  # Legend title
    ) +
    # Plot labels (clear axis and title descriptions)
    labs(
      title = title,         # Plot title (e.g., "MMP2 vs SERPINB1 expression correlation")
      x = x_label,           # X-axis label
      y = y_label,           # Y-axis label
      caption = paste("Note: Pearson correlation analysis, sample size=", nrow(data))
    ) +
    # Theme settings (clean, academic style)
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
      legend.position = "top",                               # Legend at top (avoid obscuring data)
      axis.title = element_text(face = "bold"),               # Bold axis labels
      panel.grid.major = element_line(color = "gray90"),      # Major grid lines (light gray)
      panel.grid.minor = element_blank()                      # Hide minor grid lines
    )
  
  # Add statistical results annotation (top-right corner)
  p <- p + geom_text(
    data = corr_results,
    aes(
      x = Inf,  # Right edge
      y = Inf,  # Top edge
      label = paste("r =", Pearson_r, "\nP =", P_value, "\nN =", Sample_size)  # Annotation content
    ),
    hjust = 1.1,  # Horizontal alignment (right offset)
    vjust = 1.1,  # Vertical alignment (top offset)
    size = 3.5,   # Font size
    color = "black"  # Font color
  )
  
  # Return results (plot object + correlation stats)
  return(list(plot = p, corr_stats = corr_results))
}

# 5.2 Generate correlation plots for all target genes
message("\nGenerating correlation plots...")

# Create list to store plots
cor_plots <- list()

# Generate plots for each target gene (excluding SERPINB1 itself)
for (gene in setdiff(target_genes, "SERPINB1")) {
  p <- plot_cor_scatter(
    data = merged_data,
    x_gene = "SERPINB1",
    y_gene = gene,
    title = paste(gene, "expression vs SERPINB1 methylation")
  )
  
  if (!is.null(p)) {
    cor_plots[[gene]] <- p
  }
}

# 5.3 Combine and save plots
 
  # Create grid layout - with enhanced error handling
  
    # First print individual plots to verify they render
    message("\nTesting individual plot rendering:")
    for (i in seq_along(cor_plots)) {
      print(cor_plots[[i]]$plot)
      message(sprintf(" - Plot %d/%d (%s) rendered successfully", 
                      i, length(cor_plots), names(cor_plots)[i]))
    }
    
    # Create plot grid with explicit plot objects
    plot_grid <- cowplot::plot_grid(
      plotlist = lapply(cor_plots, function(x) x$plot),  # Explicitly extract plot objects
      ncol = 2,
      labels = paste0(LETTERS[1:length(cor_plots)], ": ", names(cor_plots)),
      label_size = 12,
      align = "hv",
      axis = "l"  # Add axis alignment
    )
    
    # Add title and save with proper dimensions
    title <- ggdraw() + 
      draw_label("SERPINB1 Expression Correlations with EMT/MAPK Genes", 
                 size = 14, fontface = "bold")
    
    final_plot <- plot_grid(
      title,
      plot_grid,
      ncol = 1,
      rel_heights = c(0.1, 1)  # Title height ratio
    )
    
    # Preview before saving
    print(final_plot)
    
    # Save with dynamic sizing
    n_plots <- length(cor_plots)
    plot_rows <- ceiling(n_plots/2)
    plot_height <- max(5, 3 * plot_rows)  # Minimum 5 inches, 3 inches per row
    
    ggsave(
      filename = "SERPINB1_Expression_Correlations.pdf",
      plot = final_plot,
      width = 14,
      height = plot_height,
      dpi = 300,
      limitsize = FALSE  # Allow larger dimensions
    )
    
    # Also save as PNG for quick viewing
    ggsave(
      filename = "SERPINB1_Expression_Correlations.png",
      plot = final_plot,
      width = 14,
      height = plot_height,
      dpi = 150
    )
    
    message(sprintf(
      "\nSuccess! Saved combined plot with %d subplots to:\n- SERPINB1_Expression_Correlations.pdf\n- SERPINB1_Expression_Correlations.png",
      n_plots
    ))
    
   

# Check warnings
if (length(warnings()) > 0) {
  message("\nAdditional warnings encountered:")
  print(warnings())
}


# -------------------------- 2. save  Data --------------------------
# 2.7 Save expression matrix to CSV
write.csv(expr_matrix, "1_expression_matrix.csv", row.names = TRUE)
message("Saved raw expression matrix to: 1_expression_matrix.csv")

# 2.8 Save processed expression data
write.csv(expr_df, "2_processed_expression_data.csv", row.names = FALSE)
message("Saved processed expression data to: 2_processed_expression_data.csv")

# -------------------------- 3. Load and Preprocess Clinical Data --------------------------

# 3.5 Save clinical data
write.csv(clinical_df, "3_clinical_data_with_groups.csv", row.names = FALSE)
message("Saved clinical data with groups to: 3_clinical_data_with_groups.csv")

# -------------------------- 4. Merge Expression and Clinical Data --------------------------

# 4.2 Save merged data before NA removal
merged_data_before_na <- expr_df %>%
  inner_join(clinical_df, by = "sample_id") %>%
  mutate(across(all_of(target_genes), as.numeric))

write.csv(merged_data_before_na, "4_merged_data_before_na_removal.csv", row.names = FALSE)
message("Saved merged data before NA removal to: 4_merged_data_before_na_removal.csv")

# 4.3 Save final merged data
write.csv(merged_data, "5_final_merged_data.csv", row.names = FALSE)
message("Saved final merged data to: 5_final_merged_data.csv")

# -------------------------- 5. Correlation Analysis and Visualization --------------------------


# 5.4 Save correlation results for all genes
all_corr_results <- lapply(cor_plots, function(x) {
  x$corr_stats %>% mutate(Target_Gene = x$corr_stats$Group[1])
}) %>% bind_rows(.id = "Gene_Pair")

write.csv(all_corr_results, "6_correlation_results_all_genes.csv", row.names = FALSE)
message("Saved all correlation results to: 6_correlation_results_all_genes.csv")

# 5.5 Save individual plot data
plot_data_list <- list()
for (gene in names(cor_plots)) {
  cat(gene)
  plot_data <- merged_data %>%
    dplyr::select(sample_id, Group, SERPINB1, !!sym(gene)) %>%
    rename(Target_Gene_Expression = !!sym(gene))
  
  write.csv(plot_data, paste0("7_plot_data_SERPINB1_vs_", gene, ".csv"), row.names = FALSE)
  plot_data_list[[gene]] <- plot_data
}

for (gene in names(cor_plots)) {
  plot_data <- merged_data %>%
    dplyr::select(sample_id, Group, SERPINB1, {{gene}}) %>%
    rename(Target_Gene_Expression = {{gene}})
  
  write.csv(plot_data, paste0("7_plot_data_SERPINB1_vs_", gene, ".csv"), row.names = FALSE)
  plot_data_list[[gene]] <- plot_data
}



for (gene in names(cor_plots)) {
  plot_data <- merged_data %>%
    dplyr::select(sample_id, Group, SERPINB1, .data[[gene]]) %>%
    rename(Target_Gene_Expression = .data[[gene]])
  
  write.csv(plot_data, paste0("7_plot_data_SERPINB1_vs_", gene, ".csv"), row.names = FALSE)
  plot_data_list[[gene]] <- plot_data
}


message("Saved individual plot data to separate CSV files")

# Save combined plot data
combined_plot_data <- bind_rows(plot_data_list, .id = "Target_Gene")
write.csv(combined_plot_data, "8_all_plot_data_combined.csv", row.names = FALSE)
message("Saved combined plot data to: 8_all_plot_data_combined.csv")

# -------------------------- 6. Save Complete Analysis Metadata --------------------------
analysis_metadata <- data.frame(
  Analysis_Date = Sys.Date(),
  TCGA_Data_Version = "PRAD HiSeqV2",
  Target_Genes = paste(target_genes, collapse = ", "),
  Normal_Samples = sum(merged_data$Group == "Normal"),
  Cancer_Samples = sum(merged_data$Group == "Cancer"),
  Correlation_Method = "Pearson",
  Number_of_Plots = length(cor_plots)
)

write.csv(analysis_metadata, "9_analysis_metadata.csv", row.names = FALSE)
message("Saved analysis metadata to: 9_analysis_metadata.csv")








