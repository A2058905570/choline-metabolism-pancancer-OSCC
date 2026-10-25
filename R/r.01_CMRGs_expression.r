rm(list = ls())


setwd("./01_CMGRs_expression")
library(impute)
library(tidyverse)
library(limma)
library(GEOquery)


dat_expr <- read.csv('../00_rawdata/tcga.tpm.csv', row.names = 1)
project_sample_df <- dat_expr[, c("project", "sample_id")]
colnames(project_sample_df) <- c("cancer type", "simple_barcode")
samAnno <- project_sample_df


df <- dat_expr[, -which(names(dat_expr) == "project")]
write.csv(df, file = "dat_expr.csv", quote = TRUE, row.names = FALSE)
df <- read.csv('dat_expr.csv', row.names = 1, header = TRUE)
a <- as.data.frame(t(df))

expr_sub <- a
tumors <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
            "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
            "KIRP", "LAML", "LGG", "OV", "PCPG", "SARC",
            "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "SKCM",
            "PRAD", "READ", "STAD", "THCA", "UCEC",
            "TGCT", "THYM", "UCS", "UVM")


for (i in tumors) {
  message("--", i, "...")
  sam <- samAnno[which(samAnno$`cancer type` == i), "simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  tumsam <- comsam[substr(comsam, 14, 14) == "0"]
  norsam <- comsam[substr(comsam, 14, 14) == "1"]
  message("肿瘤样本数量: ", length(tumsam))
  message("正常样本数量: ", length(norsam))
}


rank_results <- list()

for (i in tumors) {
  message("--", i, "...")
  

  sam <- samAnno[which(samAnno$`cancer type` == i), "simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  

  tumsam <- comsam[substr(comsam, 14, 14) == "0"]
  norsam <- comsam[substr(comsam, 14, 14) == "1"]
  

  if (length(norsam) <3) {
    message("无正常样本，跳过当前肿瘤类型。")
    next
  }
  

  expr_subset <- expr_sub[, c(tumsam, norsam)]
  expr_subset[expr_subset < 0] <- 0
  expr_subset <- as.data.frame(impute.knn(as.matrix(expr_subset))$data)
  

  

  p_values <- apply(expr_subset, 1, function(x) {
    tumor_values <- x[1:length(tumsam)]
    normal_values <- x[(length(tumsam) + 1):length(x)]
    wilcox.test(tumor_values, normal_values)$p.value
  })
  

  res <- data.frame(gene = rownames(expr_subset),
                    p_value = p_values,
                    tumor = i,
                    stringsAsFactors = FALSE)
  

  res$expression <- ifelse(p_values < 0.05 & rowMeans(expr_subset[, 1:length(tumsam)]) > rowMeans(expr_subset[, (length(tumsam) + 1):length(expr_subset)]), "up",
                           ifelse(p_values < 0.05 & rowMeans(expr_subset[, 1:length(tumsam)]) < rowMeans(expr_subset[, (length(tumsam) + 1):length(expr_subset)]), "down", "no change"))
  

  res$adjusted_p_value <- p.adjust(res$p_value, method = "fdr")
  

  rank_results[[i]] <- res
}


combined_rank_results <- do.call(rbind, rank_results)


write.table(combined_rank_results, "combined_rank_sum_test_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)


head(combined_rank_results)
ndeg0<-combined_rank_results[, c("gene","tumor", "expression")]
library(dplyr)  
library(ggplot2)
library(cowplot)


df_summary <- ndeg0 %>%
  group_by(tumor, expression) %>%
  summarise(gene_count = n(), .groups = 'drop')

up_rows <- df_summary %>%
  filter(expression == "up") %>%
  mutate(row_type = "up")

down_rows <- df_summary %>%
  filter(expression == "down") %>%
  mutate(expression = "down", row_type = "down")

combined_rows <- rbind(up_rows, down_rows)

ndegs <- combined_rows[, c("tumor", "expression", "gene_count")]
colnames(ndegs) <- c("tumor", "Group", "Number")
write.csv(ndegs, file = "diff_number.csv", quote = TRUE, row.names = FALSE)


orange <- "
green <- "
library(ggplot2)
library(cowplot)
library(dplyr)


top <- ggplot(data = ndegs) +
  geom_bar(mapping = aes(x = tumor, y = Number, fill = Group),
           stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c("up" = green, "down" = orange)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1, 0, 0, 1), "lines"))
library(dplyr)

combined_rank_results$expression <- ifelse(combined_rank_results$expression %in% c("up", "down"), 
                                           combined_rank_results$expression, "other")

selected_data <- subset(combined_rank_results, 
                        !is.na(gene) & !is.na(p_value) & !is.na(expression) & !is.na(tumor) &
                          p_value < 0.05,
                        select = c(gene, p_value, expression, tumor))




colnames(selected_data) <- c("gene", "p", "expression", "tumor")


center <- ggplot(selected_data, aes(x = tumor, y = gene, fill = expression)) +
  geom_point(aes(size = -log10(p)), shape = 21, color = "black", na.rm = TRUE) +  
  scale_size_continuous(range = c(1, 5), name = "-log10(P-value)") + 
  labs(x = " ", y = " ", fill = "Expression") +
  scale_fill_manual(values = c("up" = green, "down" = orange, "other" = "#808080"), name = "Expression") +  # 设置分类颜色映射并添加图例标题
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1),  
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        legend.position = "right",  
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"))  
center

combined_plot <- plot_grid(top, center, ncol = 1, rel_heights = c(1, 4))


ggsave("01.combined_plot.pdf", plot = combined_plot, width = 10, height = 14)


ggsave("01.combined_plot.png", plot = combined_plot, width = 10, height = 14)

significant_results <- combined_rank_results[combined_rank_results$p_value < 0.05, ]
write.csv(significant_results, file = "diff_sig_pancer.csv", quote = TRUE, row.names = FALSE)


















rm(list = ls())


setwd("./01_CMGRs_expression")
library(impute)
library(tidyverse)
library(limma)
library(GEOquery)


dat_expr <- read.csv('../00_rawdata/OSCC.tpm.csv')
project_sample_df <- dat_expr[, c("project", "sample_id")]
colnames(project_sample_df) <- c("cancer type", "simple_barcode")
samAnno <- project_sample_df



df <- dat_expr[, -which(names(dat_expr) == "project")]
write.csv(df, file = "dat_OSCC.csv", quote = TRUE, row.names = FALSE)
df <- read.csv('dat_OSCC.csv', row.names = 1, header = TRUE)
a <- as.data.frame(t(df))

expr_sub <- a
tumors <- c("HNSC")


for (i in tumors) {
  message("--", i, "...")
  sam <- samAnno[which(samAnno$`cancer type` == i), "simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  tumsam <- comsam[substr(comsam, 14, 14) == "0"]
  norsam <- comsam[substr(comsam, 14, 14) == "1"]
  message("肿瘤样本数量: ", length(tumsam))
  message("正常样本数量: ", length(norsam))
}


rank_results <- list()

for (i in tumors) {
  message("--", i, "...")
  

  sam <- samAnno[which(samAnno$`cancer type` == i), "simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  

  tumsam <- comsam[substr(comsam, 14, 14) == "0"]
  norsam <- comsam[substr(comsam, 14, 14) == "1"]
  

  if (length(norsam) <3) {
    message("无正常样本，跳过当前肿瘤类型。")
    next
  }
  

  expr_subset <- expr_sub[, c(tumsam, norsam)]
  expr_subset[expr_subset < 0] <- 0
  expr_subset <- as.data.frame(impute.knn(as.matrix(expr_subset))$data)
  
  write.table(expr_subset, paste0("TCGA_", i, "_expr_subset.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
  

  p_values <- apply(expr_subset, 1, function(x) {
    tumor_values <- x[1:length(tumsam)]
    normal_values <- x[(length(tumsam) + 1):length(x)]
    wilcox.test(tumor_values, normal_values)$p.value
  })
  

  res <- data.frame(gene = rownames(expr_subset),
                    p_value = p_values,
                    tumor = i,
                    stringsAsFactors = FALSE)
  

  res$expression <- ifelse(p_values < 0.05 & rowMeans(expr_subset[, 1:length(tumsam)]) > rowMeans(expr_subset[, (length(tumsam) + 1):length(expr_subset)]), "up",
                           ifelse(p_values < 0.05 & rowMeans(expr_subset[, 1:length(tumsam)]) < rowMeans(expr_subset[, (length(tumsam) + 1):length(expr_subset)]), "down", "no change"))
  

  res$adjusted_p_value <- p.adjust(res$p_value, method = "fdr")
  

  rank_results[[i]] <- res
}


combined_rank_results <- do.call(rbind, rank_results)

significant_results <- combined_rank_results[combined_rank_results$p_value < 0.05, ]
OSCC_results <- na.omit(significant_results)
write.csv(OSCC_results, file = "diff_OSCC.csv", quote = TRUE, row.names = FALSE)
