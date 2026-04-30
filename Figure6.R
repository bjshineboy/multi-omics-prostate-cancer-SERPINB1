# 1.  ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)


# 2.  ----------------------------------------------------------
#  
data <- read.csv("/work/singleCell/prostate_singecell_20250627/GSE_Raw/sci/SERPNB1_sci_20250829/code/2/clin_methylation_data.xls", stringsAsFactors = FALSE, sep="\t")

setwd("/work/singleCell/prostate_singecell_20250627/GSE_Raw/sci/SERPNB1_sci_20250829/code/2/Figure6_Data/")
# 
meth_loci <- grep("^chr", colnames(data), value = TRUE)

# 
clinical_data <- data %>%
  dplyr::select(sampleNo, all_of(meth_loci), Urine.Collection.Timing, Age) %>%
  mutate(Urine.Collection.Timing = factor(Urine.Collection.Timing))

#  
data_long <- clinical_data %>%
  pivot_longer(
    cols = all_of(meth_loci),
    names_to = "Locus",
    values_to = "Methylation"
  ) %>%
  mutate(
    Locus = factor(Locus, levels = meth_loci),
    Urine.Collection.Timing = factor(Urine.Collection.Timing)
  )

# 3.  ----------------------------------------------------------------
timing_col <- c("Pre-operative" = "#fdae61", "Post-operative" = "#2c7bb6")

# 4.  -----------------------------------------------
#  
boxplot_list <- list()

 
for(locus in meth_loci) {
 
  locus_data <- data_long %>% filter(Locus == locus)
 
  y_max <- max(locus_data$Methylation, na.rm = TRUE) * 1.1
  
  # 
  p <- ggplot(locus_data, 
              aes(x = Urine.Collection.Timing, 
                  y = Methylation,
                  fill = Urine.Collection.Timing)) +
    geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "gray30") +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      label.y = y_max,
      size = 4,
      vjust = 1.5
    ) +
    scale_fill_manual(values = timing_col) +
    labs(
      title = paste("Methylation at", locus),
      x = "Urine Collection Timing",
      y = "Methylation Level (Beta Value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 10, angle = 0)
    )
  
  cat(locus)
  boxplot_list[[locus]] <- p
}


# 5.1  
age_boxplot_basic <- ggplot(clinical_data, 
                            aes(x = Urine.Collection.Timing, 
                                y = Age, 
                                fill = Urine.Collection.Timing)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "gray30") +
  stat_compare_means(method = "t.test", 
                     label = "p.format",
                     label.y = max(clinical_data$Age) * 1.05) +
  scale_fill_manual(values = timing_col) +
  labs(title = "Age Distribution by Urine Collection Timing",
       x = "Urine Collection Timing",
       y = "Age (Years)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none")


boxplot_list[["Age"]] <- age_boxplot_basic
# 5. 
pdf("all_loci_boxplots3.pdf", width = 12, height = 8)
for(i in seq_along(boxplot_list)) {
  print(boxplot_list[[i]])
}
dev.off()


#  
library(gridExtra)
combined_plots <- marrangeGrob(boxplot_list, nrow = 3, ncol = 3, top = NULL)
#ggsave("combined_loci_boxplots.pdf", combined_plots, width = 15, height = 10)
ggsave("combined_loci_boxplots.jpg", combined_plots, width = 15, height = 10)
# 6. 
timing_stats <- data_long %>%
  group_by(Locus, Urine.Collection.Timing) %>%
  summarise(
    Mean = mean(Methylation, na.rm = TRUE),
    Median = median(Methylation, na.rm = TRUE),
    SD = sd(Methylation, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Urine.Collection.Timing,
    values_from = c(Mean, Median, SD, N),
    names_sep = "_"   
  ) %>%
  mutate(
    Mean_Diff = `Mean_Post-operative` - `Mean_Pre-operative`,
    P_value = sapply(Locus, function(x) {
      t.test(Methylation ~ Urine.Collection.Timing, 
             data = filter(data_long, Locus == x))$p.value
    }),
    Significance = ifelse(P_value < 0.05, "*", "ns")
  )


#  
write.csv(timing_stats, "methylation_timing_statistics.csv", row.names = FALSE)

 
