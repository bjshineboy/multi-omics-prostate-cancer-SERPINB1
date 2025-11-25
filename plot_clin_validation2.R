# 1. 加载必要包 ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)

# 2. 数据准备和预处理 ----------------------------------------------------------
# 读取数据
data <- read.csv("data/clin_methylation_data.xls", stringsAsFactors = FALSE, sep="\t")

setwd("clin_validation/")
# 提取甲基化位点列名
meth_loci <- grep("^chr", colnames(data), value = TRUE)

# 准备临床数据
clinical_data <- data %>%
  select(sampleNo, all_of(meth_loci), Urine.Collection.Timing, Age) %>%
  mutate(Urine.Collection.Timing = factor(Urine.Collection.Timing))

# 转换为长格式用于箱线图
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

# 3. 颜色设置 ----------------------------------------------------------------
timing_col <- c("Pre-operative" = "#fdae61", "Post-operative" = "#2c7bb6")

# 4. 为每个位点创建分组箱线图 ------------------------------------------------
# 创建存储图形的列表
boxplot_list <- list()

# 循环生成每个位点的箱线图
for(locus in meth_loci) {
  
  # 筛选当前位点的数据
  locus_data <- data_long %>% filter(Locus == locus)
  
  # 计算y轴最大值用于调整统计标签位置
  y_max <- max(locus_data$Methylation, na.rm = TRUE) * 1.1
  
  # 创建箱线图
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


# 5.1 基础年龄分布箱线图
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
# 5. 排列并保存所有位点的箱线图 --------------------------------------------
# 方法1：将所有位点排列在一个PDF中
pdf("all_loci_boxplots3.pdf", width = 12, height = 8)
for(i in seq_along(boxplot_list)) {
  print(boxplot_list[[i]])
}
dev.off()


# 方法2：将6个位点排列在一页（2行×3列）
library(gridExtra)
combined_plots <- marrangeGrob(boxplot_list, nrow = 3, ncol = 3, top = NULL)
ggsave("combined_loci_boxplots.pdf", combined_plots, width = 15, height = 10)

# 6. 创建汇总统计表 ---------------------------------------------------------
# 计算每个位点的统计差异
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
    names_sep = "_"  # 明确指定分隔符
  ) %>%
  mutate(
    Mean_Diff = `Mean_Post-operative` - `Mean_Pre-operative`,
    P_value = sapply(Locus, function(x) {
      t.test(Methylation ~ Urine.Collection.Timing, 
             data = filter(data_long, Locus == x))$p.value
    }),
    Significance = ifelse(P_value < 0.05, "*", "ns")
  )


# 保存统计结果
write.csv(timing_stats, "methylation_timing_statistics.csv", row.names = FALSE)

# 7. 输出结果汇总 -----------------------------------------------------------
cat("分析完成！生成的文件：\n")
cat("📊 可视化结果：\n")
cat("- all_loci_boxplots.pdf：每个位点的独立箱线图\n")
cat("- combined_loci_boxplots.pdf：组合排列的箱线图（2×3布局）\n")
cat("\n📋 数据结果：\n")
cat("- methylation_timing_statistics.csv：详细的统计比较结果\n")

