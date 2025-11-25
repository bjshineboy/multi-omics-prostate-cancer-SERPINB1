# Load necessary R packages
library(rtracklayer)   # For importing bigWig files
library(GenomicRanges) # For handling genomic ranges
library(biomaRt)       # For obtaining gene coordinates
library(tidyverse)     # For data manipulation and visualization
library(patchwork)     # For combining plots
library(ggpubr)        # For adding statistical tests
library(ComplexHeatmap) # For high-quality heatmaps
library(circlize)      # For heatmap color mapping
library(showtext)      # For Chinese font support
library(colorspace)

# Set working directory
setwd("data/GSE289468_metastatic/")

# ------------------------- 1. Initialize Chinese font support -------------------------
# Load Chinese fonts
#font_add("heiti", "simhei.ttf")  # HeiTi font
#font_add("kaiti", "simkai.ttf")  # KaiTi font
showtext_auto()

# ------------------------- 2. Data preparation and preprocessing -------------------------
# Step 1: Get SERPINB1 genomic coordinates
cat("Step 1: Obtaining SERPINB1 genomic coordinates...\n")

# Use GRCh37 version (older genome assembly)
ensembl <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "grch37.ensembl.org"
)

# Get SERPINB1 genomic coordinates
serpinb1_gene <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position",
                 "transcription_start_site", "transcript_start", "transcript_end"),
  filters = "hgnc_symbol",
  values = "SERPINB1",
  mart = ensembl
)

# Create gene region GRanges objects
serpinb1_regions <- list(
  "Promoter" = GRanges(
    seqnames = paste0("chr", serpinb1_gene$chromosome_name[1]),
    ranges = IRanges(start = serpinb1_gene$transcription_start_site[1] - 2000, 
                     end = serpinb1_gene$transcription_start_site[1] + 500)
  ),
  "Gene_Body" = GRanges(
    seqnames = paste0("chr", serpinb1_gene$chromosome_name[1]),
    ranges = IRanges(start = serpinb1_gene$start_position[1], 
                     end = serpinb1_gene$end_position[1])
  ),
  "Downstream" = GRanges(
    seqnames = paste0("chr", serpinb1_gene$chromosome_name[1]),
    ranges = IRanges(start = serpinb1_gene$transcript_end[1], 
                     end = serpinb1_gene$transcript_end[1] + 2000)
  )
)

# Step 2: Prepare metadata
metadata <- data.frame(
  Sample_title = c(
    "WCM0_2, H3K27ac, ChIP-seq", "WCM0_2, H3K27me3, ChIP-seq",
    "WCM0_3, H3K27ac, ChIP-seq", "WCM0_3, H3K27me3, ChIP-seq",
    "WCM0_4, H3K27ac, ChIP-seq", "WCM0_4, H3K27me3, ChIP-seq",
    "WCM12_Z13, H3K27ac, ChIP-seq", "WCM12_Z13, H3K27me3, ChIP-seq",
    "WCM159_1, H3K27ac, ChIP-seq", "WCM159_1, H3K27me3, ChIP-seq",
    "WCM159_2, H3K27ac, ChIP-seq", "WCM159_2, H3K27me3, ChIP-seq",
    "WCM159_Z6, H3K27ac, ChIP-seq", "WCM159_Z6, H3K27me3, ChIP-seq",
    "WCM677_1, H3K27ac, ChIP-seq", "WCM677_1, H3K27me3, ChIP-seq",
    "WCM677_3, H3K27ac, ChIP-seq", "WCM677_3, H3K27me3, ChIP-seq",
    "WCM90_Z21, H3K27ac, ChIP-seq", "WCM90_Z21, H3K27me3, ChIP-seq"
  ),
  Sample_geo_accession = c(
    "GSM8791988", "GSM8791989", "GSM8791990", "GSM8791991",
    "GSM8791992", "GSM8791993", "GSM8791994", "GSM8791995",
    "GSM8791996", "GSM8791997", "GSM8791998", "GSM8791999",
    "GSM8792000", "GSM8792001", "GSM8792002", "GSM8792003",
    "GSM8792004", "GSM8792005", "GSM8792006", "GSM8792007"
  ),
  Sample_source_name_ch1 = c(
    "LN", "LN", "LN", "LN", "Pelvic mass", "Pelvic mass",
    "Brain", "Brain", "Adrenal gland", "Adrenal gland",
    "Lung", "Lung", "LN", "LN", "Bone", "Bone",
    "Liver", "Liver", "Liver", "Liver"
  ),
  stringsAsFactors = FALSE
)

# Associate bw file paths
bw_files <- paste0("./", c(
  "GSM8791988_WCM0_2_H3K27ac.bw", "GSM8791989_WCM0_2_H3K27me3.bw",
  "GSM8791990_WCM0_3_H3K27ac.bw", "GSM8791991_WCM0_3_H3K27me3.bw",
  "GSM8791992_WCM0_4_H3K27ac.bw", "GSM8791993_WCM0_4_H3K27me3.bw",
  "GSM8791994_WCM12_Z13_H3K27ac.bw", "GSM8791995_WCM12_Z13_H3K27me3.bw",
  "GSM8791996_WCM159_1_H3K27ac.bw", "GSM8791997_WCM159_1_H3K27me3.bw",
  "GSM8791998_WCM159_2_H3K27ac.bw", "GSM8791999_WCM159_2_H3K27me3.bw",
  "GSM8792000_WCM159_Z6_H3K27ac.bw", "GSM8792001_WCM159_Z6_H3K27me3.bw",
  "GSM8792002_WCM677_1_H3K27ac.bw", "GSM8792003_WCM677_1_H3K27me3.bw",
  "GSM8792004_WCM677_3_H3K27ac.bw", "GSM8792005_WCM677_3_H3K27me3.bw",
  "GSM8792006_WCM90_Z21_H3K27ac.bw", "GSM8792007_WCM90_Z21_H3K27me3.bw"
))

metadata <- metadata %>%
  mutate(
    GSM = Sample_geo_accession,
    bw_file = bw_files,
    histone_mark = str_extract(Sample_title, "H3K27ac|H3K27me3"),
    metastasis_site = case_when(
      Sample_source_name_ch1 == "LN" ~ "Lymph Node",
      Sample_source_name_ch1 == "Pelvic mass" ~ "Pelvic",
      Sample_source_name_ch1 == "Brain" ~ "Brain",
      Sample_source_name_ch1 == "Adrenal gland" ~ "Adrenal",
      Sample_source_name_ch1 == "Lung" ~ "Lung",
      Sample_source_name_ch1 == "Bone" ~ "Bone",
      Sample_source_name_ch1 == "Liver" ~ "Liver",
      TRUE ~ "Unknown"
    ),
    sample_id = str_extract(Sample_title, "^[^,]+")
  )

# ------------------------- 3. Signal extraction and analysis -------------------------
# Define signal extraction function
extract_region_signal <- function(bw_file, region_gr) {
  if (!file.exists(bw_file)) return(NA)
  tryCatch({
    bw <- import.bw(bw_file)
    overlap <- subsetByOverlaps(bw, region_gr)
    if (length(overlap) == 0) return(NA)
    mean(overlap$score, na.rm = TRUE)
  }, error = function(e) NA)
}

# Extract signals for each region
region_signals <- lapply(serpinb1_regions, function(region) {
  sapply(metadata$bw_file, function(bw_file) {
    extract_region_signal(bw_file, region)
  })
})

# Merge results
metadata <- metadata %>%
  mutate(
    promoter_signal = region_signals$Promoter,
    gene_body_signal = region_signals$Gene_Body,
    downstream_signal = region_signals$Downstream
  ) %>%
  filter(!is.na(promoter_signal) & !is.na(gene_body_signal))

# ------------------------- 4. Data visualization -------------------------
# Set theme
theme_set(theme_minimal(base_size = 12, base_family = "heiti") +
            theme(legend.position = "top",
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, color = "gray40")))

# Color scheme
plot_colors <- list(
  marks = c("H3K27ac" = "#E64B35", "H3K27me3" = "#4DBBD5"),
  sites = c("Lymph Node" = "#1B9E77", "Liver" = "#D95F02", "Bone" = "#7570B3",
            "Brain" = "#E7298A", "Lung" = "#66A61E", "Adrenal" = "#E6AB02",
            "Pelvic" = "#A6761D")
)

# A: Regional signal distribution
region_long <- metadata %>%
  pivot_longer(cols = c(promoter_signal, gene_body_signal, downstream_signal),
               names_to = "region", values_to = "signal") %>%
  mutate(region = factor(case_when(
    region == "promoter_signal" ~ "Promoter",
    region == "gene_body_signal" ~ "Gene Body",
    region == "downstream_signal" ~ "Downstream"
  ), levels = c("Promoter", "Gene Body", "Downstream")))

p1 <- ggplot(region_long, aes(x = region, y = signal, fill = histone_mark)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             aes(color = histone_mark), size = 2) +
  stat_compare_means(aes(group = histone_mark), 
                     method = "t.test", label = "p.signif", size = 5) +
  scale_fill_manual(values = plot_colors$marks) +
  scale_color_manual(values = darken(plot_colors$marks, 0.2)) +
  labs(title = "A: Epigenetic signal distribution across SERPINB1 regions",
       x = "Gene Region", y = "ChIP-seq Signal Intensity", 
       fill = "Histone Mark", color = "Histone Mark")

# B: Metastasis site-specific signals
p2 <- ggplot(metadata, aes(x = metastasis_site, y = promoter_signal, fill = histone_mark)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", alpha = 0.8) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", 
                position = position_dodge(width = 0.9), width = 0.3) +
  scale_fill_manual(values = plot_colors$marks) +
  labs(title = "B: Promoter signals by metastasis site",
       x = "Metastasis Site", y = "Promoter Signal Intensity", fill = "Histone Mark")

# C: Heatmap of regional signal patterns
heatmap_data <- metadata %>%
  dplyr::select(sample_id, metastasis_site, histone_mark, 
                promoter_signal, gene_body_signal, downstream_signal) %>%
  group_by(metastasis_site, histone_mark) %>%
  summarise(across(ends_with("signal"), mean), .groups = "drop") %>%
  pivot_longer(cols = ends_with("signal"), names_to = "region", values_to = "signal") %>%
  mutate(region = str_remove(region, "_signal")) %>%
  pivot_wider(names_from = region, values_from = signal)

mat <- as.matrix(dplyr::select(heatmap_data, promoter, gene_body, downstream))
rownames(mat) <- paste(heatmap_data$metastasis_site, heatmap_data$histone_mark, sep = "_")

ht <- Heatmap(mat, name = "Signal Intensity",
              col = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = heatmap_data$histone_mark,
              column_names_gp = gpar(fontfamily = "heiti"),
              row_names_gp = gpar(fontfamily = "heiti")
              )  

options(bitmapType = "cairo")  # Use Cairo graphics device
#dev.off()  # Close current device first
# Combine plots
final_plot <- (p1 + p2) / wrap_elements(grid.grabExpr(draw(ht))) +
  plot_annotation(title = "Epigenetic Features of SERPINB1 in Prostate Cancer Metastasis",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))


setwd("Figure5_Data")

# ---------------------------
# 1. Save Initial Metadata
# ---------------------------
write.csv(metadata, "1_initial_metadata.csv", row.names = FALSE)
cat("Initial metadata saved as: 1_initial_metadata.csv\n")

# ---------------------------
# 2. Save Gene Coordinates
# ---------------------------
write.csv(serpinb1_gene, "2_serpinb1_gene_coordinates.csv", row.names = FALSE)
cat("SERPINB1 gene coordinates saved as: 2_serpinb1_gene_coordinates.csv\n")

# ---------------------------
# 3. Save Region Definitions
# ---------------------------
region_defs <- data.frame(
  Region = names(serpinb1_regions),
  Chromosome = as.character(seqnames(serpinb1_regions[[1]])),
  Start = sapply(serpinb1_regions, function(x) start(x)),
  End = sapply(serpinb1_regions, function(x) end(x))
)
write.csv(region_defs, "3_serpinb1_region_definitions.csv", row.names = FALSE)
cat("Region definitions saved as: 3_serpinb1_region_definitions.csv\n")

# ---------------------------
# 4. Save Raw Region Signals
# ---------------------------
region_signals_df <- as.data.frame(do.call(cbind, region_signals))
region_signals_df$Sample <- rownames(region_signals_df)
write.csv(region_signals_df, "4_raw_region_signals.csv", row.names = FALSE)
cat("Raw region signals saved as: 4_raw_region_signals.csv\n")

# ---------------------------
# 5. Save Processed Metadata with Signals
# ---------------------------
write.csv(metadata, "5_processed_metadata_with_signals.csv", row.names = FALSE)
cat("Processed metadata with signals saved as: 5_processed_metadata_with_signals.csv\n")

# ---------------------------
# 6. Save Region Long Format Data
# ---------------------------
write.csv(region_long, "6_region_long_format_data.csv", row.names = FALSE)
cat("Region long format data saved as: 6_region_long_format_data.csv\n")


# ---------------------------
# 7. Save Heatmap Data
# ---------------------------
write.csv(heatmap_data, "7_heatmap_input_data.csv", row.names = FALSE)
cat("Heatmap input data saved as: 7_heatmap_input_data.csv\n")

# ---------------------------
# 8. Save Final Plot Data
# ---------------------------
plot_data <- list(
  p1_data = ggplot_build(p1)$data,
  p2_data = ggplot_build(p2)$data,
  heatmap_matrix = mat
)
saveRDS(plot_data, "8_final_plot_data.rds")
cat("Final plot data saved as: 8_final_plot_data.rds (use readRDS() to load)\n")

# ---------------------------
# 9. Save All Results in One File
# ---------------------------
all_results <- list(
  metadata = metadata,
  gene_coords = serpinb1_gene,
  region_defs = region_defs,
  region_signals = region_signals_df,
  processed_data = metadata,
  long_format = region_long,
  heatmap_data = heatmap_data
)

saveRDS(all_results, "9_all_intermediate_results.rds")
cat("All intermediate results saved as: 9_all_intermediate_results.rds\n")
cat("Use readRDS('9_all_intermediate_results.rds') to load all data\n")

# ---------------------------
# Implementation Notes:
# ---------------------------
# 1. Each CSV file is numbered for easy sequencing of the analysis pipeline
# 2. Key objects are saved at each major processing step
# 3. Both CSV (for quick viewing) and RDS (for preserving objects) formats are used
# 4. The final comprehensive save contains all intermediate results
# 5. Descriptive filenames indicate content and processing stage
# 6. Console messages confirm each save operation

# ---------------------------
# To load the saved data later:
# ---------------------------
# metadata <- read.csv("1_initial_metadata.csv")
# gene_coords <- read.csv("2_serpinb1_gene_coordinates.csv")
# all_data <- readRDS("10_all_intermediate_results.rds")


# Save plots
ggsave("SERPINB1_Epigenetic_Analysis.png", final_plot, width = 12, height = 10, dpi = 300)
ggsave("SERPINB1_Epigenetic_Analysis.pdf", final_plot, width = 12, height = 10)

cat("Analysis complete! Plots saved as SERPINB1_Epigenetic_Analysis.png and .pdf\n")


 