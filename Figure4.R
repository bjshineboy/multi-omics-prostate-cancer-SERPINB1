###############################################################################
# INTEGRATED SPATIAL TRANSCRIPTOMICS ANALYSIS OF SERPINB1 IN PROSTATE CANCER
# 
# This script performs:
# 1. Spatial analysis of SERPINB1 expression patterns
# 2. Identification of PMN-MDSC enriched regions
# 3. NE activity marker expression profiling
# 4. Integrated visualization of spatial and molecular patterns
#
# Input: 10X Visium spatial transcriptomics data
# Output: Comprehensive spatial maps and statistical analyses
###############################################################################

# ---------------------------
# 1. LOAD REQUIRED PACKAGES
# ---------------------------
library(Seurat)        # Spatial transcriptomics analysis
library(ggplot2)       # Visualization
library(dplyr)         # Data manipulation
library(ape)           # Moran's I spatial autocorrelation
library(spdep)         # Getis-Ord Gi* clustering
library(spatstat)      # Spatial point pattern analysis
library(igraph)        # Network analysis
library(viridis)       # Color schemes
library(patchwork)     # Plot arrangement
library(ggpubr)        # Publication-ready plots
library(tidyr)         # Data reshaping

# ---------------------------
# 2. DATA LOADING AND PREPROCESSING
# ---------------------------

setwd("Figure4_Data")
# Set working directory and load data
data_dir <- "data/GSE230282/"
seurat_obj <- Load10X_Spatial(data.dir = data_dir)

# Quality control - filter low-quality spots
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_Spatial > 200 & percent.mt < 20)

# Normalize data and extract expression matrix
seurat_obj <- NormalizeData(seurat_obj)
norm_expr <- GetAssayData(seurat_obj, assay = "Spatial", slot = "data")

# Extract and clean spatial coordinates
coordinates <- GetTissueCoordinates(seurat_obj) %>%
  `colnames<-`(c("x", "y", "cell_id")) %>%
  filter(complete.cases(.))

# ---------------------------
# 3. DEFINE KEY BIOMARKERS
# ---------------------------

# Target gene of interest
target_gene <- "SERPINB1"

# PMN-MDSC markers (human-specific, ranked by importance)
pmn_mdsc_markers <- c(
  "ITGAM",    # CD11b (myeloid marker)
  "FUT4",     # CD15 (granulocyte marker)
  "CEACAM8",  # CD66b (neutrophil marker)
  "CD177",    # Neutrophil activation
  "CSF3R",    # Immature granulocytes
  "S100A8",   # Immunosuppressive function
  "S100A9"    # Immunosuppressive function
)

# NE activity markers
ne_markers <- c("CXCL1", "CXCL2", "CXCL8")

# Filter available markers
available_mdsc <- intersect(pmn_mdsc_markers, rownames(norm_expr))
if(length(available_mdsc) < 1) stop("No PMN-MDSC markers found")

# ---------------------------
# 4. SPATIAL PATTERN ANALYSIS
# ---------------------------

# Calculate SERPINB1 expression thresholds
serpinb1_expr <- norm_expr[target_gene, rownames(coordinates)]
serpinb1_low_threshold <- median(serpinb1_expr, na.rm = TRUE)
serpinb1_low_cells <- names(which(serpinb1_expr < serpinb1_low_threshold))

# Calculate PMN-MDSC scores (weighted average)
mdsc_weights <- c(1, 1, 1, 0.8, 0.8, 0.8, 0.8)[1:length(available_mdsc)]
mdsc_expr <- apply(norm_expr[available_mdsc, ], 2, 
                   weighted.mean, w = mdsc_weights, na.rm = TRUE)
mdsc_high_cells <- names(which(mdsc_expr > quantile(mdsc_expr, 0.75)))

# Extract coordinates for key regions
serpinb1_low_coords <- coordinates[serpinb1_low_cells, ]
mdsc_high_coords <- coordinates[mdsc_high_cells, ]
overlap_cells <- intersect(serpinb1_low_cells, mdsc_high_cells)

# ---------------------------
# 5. SPATIAL STATISTICAL TESTS
# ---------------------------

# Moran's I spatial autocorrelation
coords_sp <- SpatialPoints(as.matrix(coordinates[, c("x", "y")]))
k <- min(10, max(4, floor(sqrt(nrow(coordinates)))))
nb <- knn2nb(knearneigh(coords_sp, k = k))
spatial_weights <- nb2listw(nb, style = "W")
moran_result <- moran.test(serpinb1_expr, spatial_weights)

# Getis-Ord Gi* local clustering
gi_star_values <- localG(serpinb1_expr, spatial_weights)

# Nearest neighbor distance analysis
if(length(serpinb1_low_cells) > 0 & length(mdsc_high_cells) > 0){
  low_ppp <- ppp(serpinb1_low_coords$x, serpinb1_low_coords$y, 
                 window = owin(range(coordinates$x), range(coordinates$y)))
  high_ppp <- ppp(mdsc_high_coords$x, mdsc_high_coords$y,
                  window = owin(range(coordinates$x), range(coordinates$y)))
  mean_dist <- mean(nncross(low_ppp, high_ppp)$dist)
  
  # Permutation test
  null_dists <- replicate(1000, {
    rand_ppp <- rpoint(npoints(high_ppp), win = as.owin(high_ppp))
    mean(nncross(low_ppp, rand_ppp)$dist)
  })
  p_val <- mean(null_dists < mean_dist)
} else {
  mean_dist <- NA
  p_val <- NA
}

# ---------------------------
# 6. DATA VISUALIZATION
# ---------------------------

# A. Create spatial distribution plot
spatial_base <- ggplot(coordinates, aes(x, y)) +
  geom_point(color = "gray90", size = 1) +
  coord_fixed() +
  theme_minimal()

spatial_plot <- spatial_base +
  # SERPINB1 low-expression regions
  geom_point(data = serpinb1_low_coords, 
             color = "#2c7bb6", size = 1.8, shape = 1, stroke = 1) +
  # PMN-MDSC high-expression regions
  geom_point(data = mdsc_high_coords, 
             color = "#d7191c", size = 1.8, shape = 4, stroke = 1) +
  # Overlapping regions
  geom_point(data = coordinates[overlap_cells, ],
             color = "black", size = 2, shape = 15) +
  labs(title = "Spatial Distribution Patterns",
       subtitle = "Blue circles: SERPINB1 low | Red crosses: PMN-MDSC high | Black squares: Overlap")

# B. Prepare NE marker expression data
ne_data <- norm_expr[ne_markers, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  left_join(coordinates, by = "cell_id") %>%
  pivot_longer(cols = all_of(ne_markers), 
               names_to = "gene", 
               values_to = "expression") %>%
  mutate(
    group = case_when(
      cell_id %in% overlap_cells ~ "Overlap",
      cell_id %in% serpinb1_low_cells ~ "SERPINB1 low",
      cell_id %in% mdsc_high_cells ~ "PMN-MDSC high",
      TRUE ~ "Other"
    ),
    group = factor(group, levels = c("SERPINB1 low", "Overlap", "PMN-MDSC high", "Other"))
  )

# C. Create violin plots for NE markers
violin_plots <- ggplot(ne_data, aes(group, expression, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  #geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = c("#2c7bb6", "#ff7f0e", "#d7191c", "gray80")) +
  facet_wrap(~gene, scales = "free_y", ncol = 1) +
  stat_compare_means(
    comparisons = list(
      c("SERPINB1 low", "Overlap"),
      c("Overlap", "PMN-MDSC high"),
      c("PMN-MDSC high", "Other")
    ),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  labs(x = "Spatial Region", y = "Expression Level") +
  theme_minimal()





# D. Combine plots
final_plot <- (spatial_plot + labs(tag = "A")) / 
  (violin_plots + labs(tag = "B")) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title = "Integrated Spatial and Molecular Analysis of SERPINB1",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

# ---------------------------
# 7. SAVE RESULTS
# ---------------------------

# Save combined plot
ggsave("Integrated_SERPINB1_Analysis.png", 
       final_plot, 
       width = 12, height = 14, dpi = 300)

# Save statistical results
write.csv(data.frame(
  Metric = c("Moran's I", "P-value (Moran's)", "Mean NN Distance", "P-value (NN)"),
  Value = c(moran_result$estimate[1], moran_result$p.value, mean_dist, p_val)
), "spatial_statistics.csv")

# Print summary
cat("\n===== ANALYSIS COMPLETED =====\n")
cat("- Moran's I:", round(moran_result$estimate[1], 3), "\n")
cat("- Spatial co-localization p-value:", format.pval(p_val), "\n")
cat("- Results saved to:\n")
cat("  Integrated_SERPINB1_Analysis.png\n")
cat("  spatial_statistics.csv\n")







###############################################################################
# SPATIAL TRANSCRIPTOMICS ANALYSIS PIPELINE WITH INTERMEDIATE FILE SAVING
###############################################################################

# ---------------------------
# 1. INITIALIZATION AND SETUP
# ---------------------------
# Create output directory structure
output_dirs <- c("raw_data", "processed_data", "quality_control", 
                 "spatial_stats", "expression_profiles", "final_results")
sapply(output_dirs, function(x) if(!dir.exists(x)) dir.create(x))

# Initialize logging
log_file <- file("analysis_log.txt", open = "wt")
sink(log_file, type = "output", split = TRUE)
cat(format(Sys.time(), "\n===== ANALYSIS STARTED: %Y-%m-%d %H:%M:%S =====\n"))

# ---------------------------
# 2. DATA LOADING AND QC
# ---------------------------
cat("\n--- Loading and QC ---\n")
 
  # Load spatial data
  seurat_obj <- Load10X_Spatial(data.dir = data_dir)
  cat("Raw data dimensions:", dim(seurat_obj), "\n")
  
  # Save raw counts
  
  counts_matrix <- LayerData(seurat_obj, assay = "Spatial", layer = "counts")
  saveRDS(counts_matrix, 
          "raw_data/raw_counts_matrix.rds")
  write.csv(GetTissueCoordinates(seurat_obj), 
            "raw_data/raw_spatial_coordinates.csv")
  
  # QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  qc_metrics <- data.frame(
    nFeature = seurat_obj$nFeature_Spatial,
    nCount = seurat_obj$nCount_Spatial,
    percent.mt = seurat_obj$percent.mt
  )
  write.csv(qc_metrics, "quality_control/pre_filter_qc_metrics.csv")
  
  # Filtering
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_Spatial > 200 & percent.mt < 20)
  cat("After QC filtering:", dim(seurat_obj), "\n")
  write.csv(qc_metrics[colnames(seurat_obj),], 
            "quality_control/post_filter_qc_metrics.csv")
  


# ---------------------------
# 3. DATA NORMALIZATION
# ---------------------------
cat("\n--- Normalization ---\n")
seurat_obj <- NormalizeData(seurat_obj)
norm_data <- GetAssayData(seurat_obj, slot = "data")
write.csv(norm_data, "processed_data/normalized_expression_matrix.csv", row.names = TRUE)

# Save normalized expression for key genes
key_genes <- c(target_gene, pmn_mdsc_markers, ne_markers)
write.csv(norm_data[key_genes, ], 
          "processed_data/normalized_key_genes.csv")

# ---------------------------
# 4. SPATIAL COORDINATES PROCESSING
# ---------------------------
coordinates <- GetTissueCoordinates(seurat_obj) %>%
  `colnames<-`(c("x", "y", "cell_id")) %>%
  filter(complete.cases(.))
write.csv(coordinates, "processed_data/cleaned_spatial_coordinates.csv")

# ---------------------------
# 5. TARGET GENE ANALYSIS
# ---------------------------
cat("\n--- Target Gene Analysis ---\n")
serpinb1_expr <- norm_data[target_gene, rownames(coordinates)]
serpinb1_stats <- data.frame(
  cell_id = names(serpinb1_expr),
  expression = as.numeric(serpinb1_expr),
  stringsAsFactors = FALSE
) %>%
  left_join(coordinates, by = "cell_id")

write.csv(serpinb1_stats, "expression_profiles/serpinb1_expression_profiles.csv")

# ---------------------------
# 6. PMN-MDSC SCORE CALCULATION
# ---------------------------
cat("\n--- PMN-MDSC Scoring ---\n")
mdsc_weights <- c(1, 1, 1, 0.8, 0.8, 0.8, 0.8)[1:length(available_mdsc)]
mdsc_scores <- apply(norm_data[available_mdsc, ], 2, 
                     weighted.mean, w = mdsc_weights, na.rm = TRUE)

mdsc_score_df <- data.frame(
  cell_id = names(mdsc_scores),
  mdsc_score = mdsc_scores,
  stringsAsFactors = FALSE
) %>%
  left_join(coordinates, by = "cell_id")

write.csv(mdsc_score_df, "expression_profiles/pmn_mdsc_scores.csv")

# ---------------------------
# 7. SPATIAL STATISTICS
# ---------------------------
cat("\n--- Spatial Statistics ---\n")

# Moran's I
moran_result <- moran.test(serpinb1_expr, spatial_weights)
moran_df <- data.frame(
  Statistic = c("Moran's I", "Expected", "Variance", "P-value"),
  Value = c(moran_result$estimate[1], 
            moran_result$estimate[2],
            moran_result$estimate[3],
            moran_result$p.value)
)
write.csv(moran_df, "spatial_stats/morans_i_results.csv", row.names = FALSE)

# Gi* Analysis
# 先计算所有需要的向量
p_values <- 2 * pnorm(-abs(gi_star_values))
significance <- ifelse(abs(gi_star_values) > 1.96 & 
                         p_values < 0.05, 
                       ifelse(gi_star_values > 0, "Hotspot", "Coldspot"), 
                       "Not significant")

# 然后构建data.frame
gi_star_results <- data.frame(
  cell_id = names(gi_star_values),
  Gi_star = as.numeric(gi_star_values),
  p_value = p_values,
  significance = significance
)

write.csv(gi_star_results, "spatial_stats/getis_ord_gi_results.csv", row.names = FALSE)

# ---------------------------
# 8. NE ACTIVITY MARKER ANALYSIS
# ---------------------------
cat("\n--- NE Activity Analysis ---\n")
ne_expr <- norm_data[ne_markers, ]
ne_expr_long <- as.data.frame(t(ne_expr)) %>%
  tibble::rownames_to_column("cell_id") %>%
  pivot_longer(cols = -cell_id, names_to = "gene", values_to = "expression") %>%
  left_join(coordinates, by = "cell_id")

write.csv(ne_expr_long, "expression_profiles/ne_activity_markers_long.csv", row.names = FALSE)

# ---------------------------
# 9. REGION DEFINITION AND OVERLAP ANALYSIS
# ---------------------------
cat("\n--- Region Definition ---\n")

# Define regions
region_def <- data.frame(
  cell_id = rownames(coordinates),
  serpinb1_low = ifelse(rownames(coordinates) %in% serpinb1_low_cells, 1, 0),
  mdsc_high = ifelse(rownames(coordinates) %in% mdsc_high_cells, 1, 0),
  overlap = ifelse(rownames(coordinates) %in% overlap_cells, 1, 0)
)

write.csv(region_def, "processed_data/region_definitions.csv", row.names = FALSE)

# ---------------------------
# 10. INTEGRATED VISUALIZATION
# ---------------------------
cat("\n--- Visualization ---\n")

# Save individual plots
ggsave("final_results/spatial_distribution.png", spatial_plot, 
       width = 8, height = 6, dpi = 300)
ggsave("final_results/ne_marker_violins.png", violin_plots, 
       width = 8, height = 10, dpi = 300)
ggsave("final_results/combined_analysis.png", final_plot, 
       width = 12, height = 14, dpi = 300)

# ---------------------------
# 11. FINAL OUTPUT AND CLEANUP
# ---------------------------
# Create analysis report
report_data <- list(
  Parameters = data.frame(
    Parameter = c("Target Gene", "PMN-MDSC Markers", "NE Markers",
                  "Moran's I", "Gi* Hotspots", "Gi* Coldspots"),
    Value = c(target_gene, 
              paste(available_mdsc, collapse = ", "),
              paste(ne_markers, collapse = ", "),
              round(moran_result$estimate[1], 3),
              sum(gi_star_results$significance == "Hotspot"),
              sum(gi_star_results$significance == "Coldspot"))
  ),
  SummaryStats = data.frame(
    Statistic = c("Total Spots", "SERPINB1 Low", "PMN-MDSC High", "Overlap"),
    Count = c(ncol(seurat_obj),
              length(serpinb1_low_cells),
              length(mdsc_high_cells),
              length(overlap_cells))
  )
)

saveRDS(report_data, "final_results/analysis_summary.rds")

# Zip intermediate files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
zip(paste0("intermediate_files_", timestamp, ".zip"), 
    unlist(sapply(output_dirs, list.files, full.names = TRUE)))

cat("\n===== ANALYSIS COMPLETED SUCCESSFULLY =====\n")
cat("Intermediate files saved to:", paste0("intermediate_files_", timestamp, ".zip"), "\n")
cat("Final results saved to final_results/\n")

sink()
close(log_file)












