# --------------------------- 0. Load R Packages ---------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(ggpubr)

setwd("Figure3_Data")

# --------------------------- 1. Data Loading and QC ---------------------------
seurat_obj <- Read10X_h5("data/PRAD_GSE137829_expression.h5") %>% 
  CreateSeuratObject(
    project = "Prostate_Cancer",
    min.cells = 3,
    min.features = 200,
    names.field = 1,
    names.delim = "-"
  )

# QC: Calculate mitochondrial and ribosomal percentages
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Filter low-quality cells
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 6000 &
                       percent.mt < 10 &
                       percent.rb < 40)

# --------------------------- 2. Data Preprocessing ---------------------------
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("percent.mt"))

# --------------------------- 3. Dimensionality Reduction and Clustering ---------------------------
seurat_obj <- seurat_obj %>%
  RunPCA(npcs = 50) %>%
  RunTSNE(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = c(0.3, 0.5, 0.8))

# --------------------------- 4. Club-like Cell Identification ---------------------------
club_markers <- c("LTF", "SCGB3A1", "KRT5", "KRT15")
seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(Club_signature = club_markers),
                             name = "Club_score")

club_threshold <- quantile(seurat_obj$Club_score1, 0.8)
seurat_obj$cell_type <- ifelse(seurat_obj$Club_score1 > club_threshold, 
                               "Club-like", "Other")
Idents(seurat_obj) <- "cell_type"

# --------------------------- 5. t-SNE Visualization ---------------------------
cell_types <- unique(seurat_obj$cell_type)
color_list <- ifelse(cell_types == "Club-like", "red", "gray85")
names(color_list) <- cell_types

tsne_plot_club <- DimPlot(seurat_obj, 
                          reduction = "tsne", 
                          group.by = "cell_type", 
                          label = FALSE, 
                          pt.size = 0.4) +
  scale_color_manual(values = color_list) +
  ggtitle("t-SNE Plot: Club-like Cells (Red)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom",
        legend.title = element_blank())

# Calculate Club-like percentage
club_count <- sum(seurat_obj$cell_type == "Club-like")
total_cells <- nrow(seurat_obj@meta.data)
club_percent <- round((club_count / total_cells) * 100, 1)

tsne_final <- tsne_plot_club +
  annotate("text", 
           x = max(seurat_obj@reductions$tsne@cell.embeddings[,1]) * 0.8, 
           y = max(seurat_obj@reductions$tsne@cell.embeddings[,2]) * 0.9,
           label = paste0(club_percent, "% of Tumor Cells"),
           color = "red", 
           size = 5, 
           fontface = "bold")


# ---------------------------- SERPINB1 Mrna in club-like and other
# Extract SERPINB1 expression data with cell type information
serpinb1_data <- data.frame(
  Expression = GetAssayData(seurat_obj, slot = "data")["SERPINB1", ],
  CellType = seurat_obj$cell_type
) %>%
  # Remove any NA values
  na.omit() %>%
  filter(Expression !=0.0 ) %>%
  # Convert to factors for proper ordering
  mutate(CellType = factor(CellType, levels = c("Club-like", "Other")))

# Calculate summary statistics for annotation
summary_stats <- serpinb1_data %>%
  group_by(CellType) %>%
  summarise(
    Median = median(Expression),
    Q1 = quantile(Expression, 0.25),
    Q3 = quantile(Expression, 0.75),
    .groups = 'drop'
  )

# Create boxplot with jitter points
serpinb1_boxplot <- ggplot(serpinb1_data, aes(x = CellType, y = Expression)) +
  # Boxplot layer
  geom_boxplot(
    aes(fill = CellType),
    width = 0.3,
    outlier.shape = NA, # Hide default outliers (we'll show all points)
    alpha = 0.7
  )  +
  # Statistical comparison
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Club-like", "Other")),
    label = "p.format",
    label.y = max(serpinb1_data$Expression) * 1.1,
    size = 4
  ) +
  # Custom colors
  scale_fill_manual(values = c("Club-like" = "red", "Other" = "gray85")) +
  scale_color_manual(values = c("Club-like" = "red", "Other" = "gray85")) +
  # Labels and theme
  labs(
    title = "SERPINB1 Expression by Cell Type",
    x = "",
    y = "Log-normalized Expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 10)
  ) +
  # Add sample size annotation
  annotate(
    "text",
    x = 1:2,
    y = min(serpinb1_data$Expression) * 0.95,
    label = paste0("n=", table(serpinb1_data$CellType)),
    size = 3.5
  )

# Save the plot
ggsave("SERPINB1_boxplot_jitter.pdf", 
       serpinb1_boxplot,
       width = 5, 
       height = 6,
       dpi = 300)

# Show the plot
print(serpinb1_boxplot)


# --------------------------- 6. NE Activity Analysis ---------------------------
ne_markers <- c("CXCL1", "CXCL2", "CXCL8")
seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(NE_activity = ne_markers),
                             name = "NE_score")

# --------------------------- 6.1 Compare NE Activity Between Club-like and Other Cells ---------------------------
# Create boxplot with jitter points for NE activity
ne_boxplot <- ggplot(seurat_obj@meta.data, 
                     aes(x = cell_type, y = NE_score1)) +
  # Boxplot layer
  geom_boxplot(
    aes(fill = cell_type),
    width = 0.3,
    outlier.shape = NA,  # We'll show all points via jitter
    alpha = 0.7
  ) +

  # Statistical comparison
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Club-like", "Other")),
    label = "p.format",
    label.y = max(seurat_obj$NE_score1, na.rm = TRUE) * 1.1,
    size = 4
  ) +
  # Color scheme matching t-SNE plot
  scale_fill_manual(values = color_list) +
  scale_color_manual(values = color_list) +
  # Labels and title
  labs(
    title = "NE Activity: Club-like vs Other Cells",
    y = "NE Activity Score (NE_score1)",
    x = "Cell Type"
  ) +
  # Theme adjustments
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "none"
  ) +
  # Add sample size annotation
  annotate(
    "text",
    x = 1:2,
    y = min(seurat_obj$NE_score1, na.rm = TRUE) * 0.95,
    label = paste0("n=", table(seurat_obj$cell_type)),
    size = 3.5
  )

# Calculate summary statistics for annotation
summary_ne_stats <- seurat_obj@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    Median = median(NE_score1),
    Q1 = quantile(NE_score1, 0.25),
    Q3 = quantile(NE_score1, 0.75),
    .groups = 'drop'
  )

summary_ne_stats
# Save the plot
ggsave("NE_activity_boxplot.pdf", 
       ne_boxplot,
       width = 6, 
       height = 6)

 
# --------------------------- 8. Differential Expression Analysis ---------------------------
club_markers_de <- FindMarkers(seurat_obj,
                               ident.1 = "Club-like",
                               ident.2 = "Other",
                               logfc.threshold = 0.25,
                               min.pct = 0.1,
                               only.pos = TRUE)

sig_genes <- club_markers_de %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

# --------------------------- 9. Functional Enrichment Analysis ---------------------------
ego <- enrichGO(gene = sig_genes$gene,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

kk <- enrichKEGG(gene = bitr(sig_genes$gene, 
                             fromType = "SYMBOL",
                             toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)



write.csv(ego@result, "ego_result.csv")
write.csv(kk@result, "kk_result.csv")
# Visualize enrichment results
go_plot <- dotplot(ego, showCategory = 15) + 
  ggtitle("GO Biological Process Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12))

kegg_plot <- dotplot(kk, showCategory = 15) + 
  ggtitle("KEGG Pathway Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12))

# --------------------------- 10. SERPINB1 Correlation Analysis ---------------------------
club_cells <- subset(seurat_obj, subset = cell_type == "Club-like")
target_gene <- "SERPINB1"
ne_markers <- c("CXCL1", "CXCL2", "CXCL8")

# Extract expression matrix (log-normalized)
expr_matrix <- GetAssayData(club_cells, assay = "RNA", slot = "data")

# Prepare data frame for correlation analysis
data_for_corr <- data.frame(
  SERPINB1 = as.numeric(expr_matrix[target_gene, ]),
  cell_id = colnames(club_cells)
)

# Add NE marker expression
for (marker in ne_markers) {
  marker_expr <- tryCatch(
    expr_matrix[marker, ], 
    error = function(e) rep(NA, ncol(club_cells)) # Handle missing genes
  )
  data_for_corr[[marker]] <- as.numeric(marker_expr)
}

# Remove NA values
data_for_corr <- na.omit(data_for_corr) %>%
  filter(SERPINB1 != 0.0 ) %>%
  filter(CXCL1 != 0.0 ) %>%
  filter(CXCL2 != 0.0 ) %>%
  filter(CXCL8 != 0.0 )

# Create correlation plots for each NE marker
plot_list <- list()
for (marker in ne_markers) {
  # Calculate Spearman correlation (non-parametric)
  spearman_cor <- cor(data_for_corr$SERPINB1, data_for_corr[[marker]], method = "spearman")
  spearman_pval <- cor.test(data_for_corr$SERPINB1, data_for_corr[[marker]], method = "spearman")$p.value
  
  # Generate plot
  plot_list[[marker]] <- ggplot(data_for_corr, aes(x = SERPINB1, y = .data[[marker]])) +
    geom_point(alpha = 0.3, size = 1.2, color = "#3498db") + # Blue points (translucent)
    geom_smooth(method = "lm", se = FALSE, color = "#e74c3c", linewidth = 1) + # Red trend line
    labs(
      title = paste("SERPINB1 vs", marker),
      x = "SERPINB1 Expression (log-normalized)",
      y = paste(marker, "Expression (log-normalized)")
    ) +
    annotate(
      "text", 
      x = min(data_for_corr$SERPINB1), 
      y = max(data_for_corr[[marker]]),
      label = sprintf("ρ = %.2f\nP = %.2e", spearman_cor, spearman_pval),
      hjust = 0, vjust = 1, size = 4 # Position label at top-left
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 10)
    )
}

# Combine plots into a single panel
serpinb1_corr_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "SERPINB1 vs NE Markers in Club-like Cells",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

# --------------------------- 11. Combine All Plots Into Final Figure ---------------------------
final_plot <- (tsne_final | ne_boxplot | serpinb1_boxplot) / # Row 1: t-SNE + Violin plot (Club-like features)
  #p_combined / # Row 2: CXCL vs NE activity (NE functional output)
  (go_plot | kegg_plot) / # Row 3: GO + KEGG (Differential gene function)
  serpinb1_corr_plot + # Row 4: SERPINB1 vs NE markers (Molecular mechanism)
  plot_annotation(
    title = "Comprehensive Analysis of Club-like Cells in Prostate Cancer",
    tag_levels = "A", # Add panel tags (A, B, C, D)
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.margin = margin(10, 10, 10, 10) # Add margin around the figure
    )
  ) +
  plot_layout(
    heights = c(1, 1.2, 1.2, 1), # Adjust row heights (relative)
    widths = c(1, 1), # Set equal width for Row 1 plots
    guides = "collect" # Collect legends (if any) to avoid duplication
  )

 
# --------------------------- 12. Save Results ---------------------------


# Save final combined figure (high-resolution PDF)
ggsave(
  filename = "final_combined_figure.png",
  plot = final_plot,
  width = 18, # Increased width for Row 1 (two plots)
  height = 20, # Increased height for four rows
  dpi = 300 # High resolution for publication
)

# Save individual plots (optional)
ggsave("tsne_club.pdf", tsne_final, width = 8, height = 6, dpi = 300)
ggsave("ne_activity_boxplot.pdf", ne_boxplot, width = 8, height = 6, dpi = 300)
ggsave("cxcl_ne_correlation2.pdf", p_combined, width = 12, height = 6, dpi = 300)
ggsave("go_enrichment.pdf", go_plot, width = 8, height = 6, dpi = 300)
ggsave("kegg_enrichment.pdf", kegg_plot, width = 8, height = 6, dpi = 300)
ggsave("serpinb1_ne_correlation.pdf", serpinb1_corr_plot, width = 12, height = 6, dpi = 300)

# Save Seurat object and differential genes
saveRDS(seurat_obj, "prostate_seurat.rds")
write.csv(sig_genes, "club_cell_markers_sig_genes.csv", row.names = FALSE)

# --------------------------- 13. Results Summary ---------------------------
cat("\n=== ANALYSIS SUMMARY ===\n")
cat(sprintf("1. Club-like cells count: %d\n", club_count))
cat(sprintf("2. Club-like cells percentage: %.1f%%\n", club_percent))
cat(sprintf("3. Significant upregulated genes in Club-like cells: %d\n", nrow(sig_genes)))

ne_pvalue <- wilcox.test(NE_score1 ~ cell_type, data = seurat_obj@meta.data)$p.value
cat(sprintf("4. NE activity in Club-like cells vs other cells: p = %.3f (Wilcoxon test)\n", ne_pvalue))

cat("5. All plots saved to: \n")








########################## test ################

# --------------------------- 10. SERPINB1 Correlation Analysis ---------------------------
other_cells <- subset(seurat_obj, subset = cell_type == "Other")
target_gene <- "SERPINB1"
ne_markers <- c("CXCL1", "CXCL2", "CXCL8")

# Extract expression matrix (log-normalized)
expr_matrix_other <- GetAssayData(other_cells, assay = "RNA", slot = "data")

# Prepare data frame for correlation analysis
data_for_corr_other <- data.frame(
  SERPINB1 = as.numeric(expr_matrix_other[target_gene, ]),
  cell_id = colnames(other_cells)
)

# Add NE marker expression
for (marker in ne_markers) {
  marker_expr <- tryCatch(
    expr_matrix_other[marker, ], 
    error = function(e) rep(NA, ncol(other_cells)) # Handle missing genes
  )
  data_for_corr_other[[marker]] <- as.numeric(marker_expr)
}

# Remove NA values
data_for_corr_other <- na.omit(data_for_corr_other) %>%
  filter(SERPINB1 != 0.0 ) %>%
  filter(CXCL1 != 0.0 ) %>%
  filter(CXCL2 != 0.0 ) %>%
  filter(CXCL8 != 0.0 )

# Create correlation plots for each NE marker
plot_list2 <- list()
for (marker in ne_markers) {
  # Calculate Spearman correlation (non-parametric)
  spearman_cor <- cor(data_for_corr_other$SERPINB1, data_for_corr_other[[marker]], method = "spearman")
  spearman_pval <- cor.test(data_for_corr_other$SERPINB1, data_for_corr_other[[marker]], method = "spearman")$p.value
  
  # Generate plot
  plot_list2[[marker]] <- ggplot(data_for_corr_other, aes(x = SERPINB1, y = .data[[marker]])) +
    geom_point(alpha = 0.3, size = 1.2, color = "#3498db") + # Blue points (translucent)
    geom_smooth(method = "lm", se = FALSE, color = "#e74c3c", linewidth = 1) + # Red trend line
    labs(
      title = paste("SERPINB1 vs", marker),
      x = "SERPINB1 Expression (log-normalized)",
      y = paste(marker, "Expression (log-normalized)")
    ) +
    annotate(
      "text", 
      x = min(data_for_corr_other$SERPINB1), 
      y = max(data_for_corr_other[[marker]]),
      label = sprintf("ρ = %.2f\nP = %.2e", spearman_cor, spearman_pval),
      hjust = 0, vjust = 1, size = 4 # Position label at top-left
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 10)
    )
}

# Combine plots into a single panel
serpinb1_corr_plot_other <- wrap_plots(plot_list2, ncol = 3) +
  plot_annotation(
    title = "SERPINB1 vs NE Markers in other Cells",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )


serpinb1_corr_plot_other








