# ====================== Gene Expression Correlation Analysis (TCGA-PRAD) ======================
# Objective: Analyze expression correlations between SERPINB1 and EMT/MAPK genes
# Data sources: TCGA-PRAD mRNA expression data (HiSeqV2)
# Output: Combined scatter plot matrix
# ==============================================================================

# -------------------------- 1. Load Required Packages --------------------------
library(data.table)  # Efficient large file reading
library(dplyr)       # Data cleaning
library(tidyr)       # Data reshaping
library(ggplot2)     # Plotting
library(cowplot)     # Multi-panel layouts
library(broom)       # Tidy statistical results
library(ggrepel)     # Prevent text overlap

setwd("/work/singleCell/prostate_singecell_20250627/GSE_Raw/sci/SERPNB1_sci_20250829/code/2/Figure4_Data/")

# -------------------------- 2. Load and Preprocess Data --------------------------
# Define target genes
target_genes <- c("SERPINB1", "MMP2", "MMP9", "SNAI1", "TWIST1", "MAPK3", "MAPK1")

# Load expression data
expr_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_HiSeqV2.gz"
expr_data <- data.table::fread(
  cmd = paste("zcat", expr_file),
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

# Process expression matrix
gene_names <- expr_data[-1, 1] %>% as.character()
sample_names <- expr_data[1, -1] %>% as.character()
expr_matrix <- as.matrix(expr_data[-1, -1])
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- sample_names

# Convert to tidy format
expr_df <- as.data.frame(expr_matrix) %>%
  tibble::rownames_to_column("Gene") %>%
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

# Load clinical data
clinical_file <- "/work/singleCell/prostate_singecell_20250627/GSE_Raw/DOWN_FROM_LJY/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix"
clinical_data <- data.table::fread(
  file = clinical_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

# Process clinical data
clinical_df <- clinical_data %>%
  mutate(sample_id = clinical_data[["sampleID"]]) %>%
  dplyr::select(sample_id, all_of("sample_type")) %>%
  rename(sample_type = all_of("sample_type")) %>%
  mutate(
    Group = case_when(
      sample_type == "Solid Tissue Normal" ~ "Normal",
      grepl("Tumor|Primary", sample_type, ignore.case = TRUE) ~ "Cancer",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Group))

# Merge data
merged_data <- expr_df %>%
  inner_join(clinical_df, by = "sample_id") %>%
  mutate(across(all_of(target_genes), as.numeric)) %>%
  na.omit()

# -------------------------- 3. Define Plotting Function --------------------------
plot_cor_scatter <- function(data, x_gene, y_gene, group_var = "Group",
                             x_label = paste(x_gene, "expression"),
                             y_label = paste(y_gene, "expression"),
                             title = paste(y_gene, "vs", x_gene),
                             color_palette = c("Normal" = "#0072B2", "Cancer" = "#D55E00")) {
  
  # Calculate Pearson correlation
  corr_results <- data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      n = sum(!is.na(!!sym(x_gene)) & !is.na(!!sym(y_gene))),
      cor_test = list(cor.test(!!sym(x_gene), !!sym(y_gene), method = "pearson"))
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
      Pearson_r = round(Pearson_r, 3),
      P_value_format = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value))
    )
  
  #  
  x_range <- range(data[[x_gene]], na.rm = TRUE)
  y_range <- range(data[[y_gene]], na.rm = TRUE)
  x_pos <- x_range[1] + 0.7 * diff(x_range)  #  
  y_pos_normal <- y_range[1] + 0.2 * diff(y_range)  # 
  y_pos_cancer <- y_range[1] + 0.8 * diff(y_range)  
  
  # Create scatter plot
  p <- ggplot(data, aes(
    x = !!sym(x_gene),
    y = !!sym(y_gene),
    color = !!sym(group_var)
  )) +
    geom_point(alpha = 0.6, size = 1.5) +  #  
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      linetype = "dashed",
      size = 0.5
    ) +
    scale_color_manual(
      values = color_palette,
      name = "Group"
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),  #  
      legend.position = "top",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "pt")  # 
    )
  
  #  
  for (i in 1:nrow(corr_results)) {
    group <- corr_results$Group[i]
    y_pos <- ifelse(group == "Normal", y_pos_normal, y_pos_cancer)
    
    p <- p + annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = sprintf("%s: r = %.3f\nP = %s\nn = %d", 
                      group, 
                      corr_results$Pearson_r[i], 
                      corr_results$P_value_format[i],
                      corr_results$Sample_size[i]),
      size = 2.5,  # 
      hjust = 0,
      vjust = 0.5,
      color = "black"
    )
  }
  
  return(list(plot = p, corr_stats = corr_results))
}

# -------------------------- 4. Generate Combined Plot --------------------------
# Generate plots for each target gene (excluding SERPINB1 itself)
cor_plots <- list()
all_corr_results <- list()

for (gene in setdiff(target_genes, "SERPINB1")) {
  result <- plot_cor_scatter(
    data = merged_data,
    x_gene = "SERPINB1",
    y_gene = gene,
    title = gene
  )
  
  if (!is.null(result)) {
    cor_plots[[gene]] = result$plot
    #  
    all_corr_results[[gene]] <- result$corr_stats %>%
      mutate(Gene = gene, Comparison = paste(gene, "vs SERPINB1"))
  }
}

# 
if (length(all_corr_results) > 0) {
  corr_summary <- bind_rows(all_corr_results) %>%
    dplyr::select(Gene, Comparison, Group, Pearson_r, P_value, Sample_size) %>%
    arrange(Gene, Group)
  
  write.csv(corr_summary, "SERPINB1_Expression_Correlations_Results.csv", row.names = FALSE)
  cat("Saved correlation results to: SERPINB1_Expression_Correlations_Results.csv\n")
  
  # 
  write.csv(merged_data, "Processed_Expression_Data.csv", row.names = FALSE)
  cat("Saved processed expression data to: Processed_Expression_Data.csv\n")
}

#  
group_stats <- merged_data %>%
  group_by(Group) %>%
  summarise(
    n_samples = n(),
    SERPINB1_mean = mean(SERPINB1, na.rm = TRUE),
    SERPINB1_sd = sd(SERPINB1, na.rm = TRUE)
  )

write.csv(group_stats, "Group_Statistics.csv", row.names = FALSE)
cat("Saved group statistics to: Group_Statistics.csv\n")

# Create plot grid
plot_grid_obj <- cowplot::plot_grid(
  plotlist = cor_plots,
  ncol = 2,
  labels = paste0(LETTERS[1:length(cor_plots)], ": ", names(cor_plots)),
  label_size = 10,  
  label_fontface = "bold",
  align = "hv",
  axis = "l",
  hjust = 0,
  vjust = 1.1,  #  
  label_x = 0.05,  #  
  label_y = 0.95    
)

# Add title
title <- ggdraw() + 
  draw_label("SERPINB1 Expression Correlations with EMT/MAPK Genes", 
             size = 12, fontface = "bold")   

#  
final_plot <- plot_grid(
  title,
  plot_grid_obj,
  ncol = 1,
  rel_heights = c(0.07, 1)   
)

# 
n_plots <- length(cor_plots)
plot_rows <- ceiling(n_plots/2)
#  
plot_height <- max(8, 3.5 * plot_rows)   
plot_width <- 12   

 
ggsave(
  filename = "SERPINB1_Expression_Correlations.pdf",
  plot = final_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)
 
ggsave(
  filename = "SERPINB1_Expression_Correlations.png",
  plot = final_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

 
analysis_summary <- data.frame(
  Analysis_Date = Sys.Date(),
  Dataset = "TCGA-PRAD",
  Target_Genes = paste(target_genes, collapse = ", "),
  Normal_Samples = sum(merged_data$Group == "Normal"),
  Cancer_Samples = sum(merged_data$Group == "Cancer"),
  Total_Samples = nrow(merged_data),
  Number_of_Plots = length(cor_plots)
)

write.csv(analysis_summary, "Analysis_Summary.csv", row.names = FALSE)
cat("Saved analysis summary to: Analysis_Summary.csv\n")

 
 
cat("ANALYSIS SUMMARY\n")

cat(sprintf("Analysis Date: %s\n", Sys.Date()))
cat(sprintf("Total Samples: %d\n", nrow(merged_data)))
cat(sprintf("Normal Samples: %d\n", sum(merged_data$Group == "Normal")))
cat(sprintf("Cancer Samples: %d\n", sum(merged_data$Group == "Cancer")))
cat(sprintf("Number of Gene Comparisons: %d\n", length(cor_plots)))
 

cat("Success! Generated the following files:\n")
cat("1. SERPINB1_Expression_Correlations.pdf - Combined correlation plot\n")
cat("2. SERPINB1_Expression_Correlations.png - Image version\n")
cat("3. SERPINB1_Expression_Correlations_Results.csv - Correlation coefficients\n")
cat("4. Processed_Expression_Data.csv - Cleaned expression data\n")
cat("5. Group_Statistics.csv - Group-wise statistics\n")
cat("6. Analysis_Summary.csv - Analysis metadata\n")


