rm(list = ls())


setwd("./07_pan_GSVA")
gene<-read.csv('../00_rawdata/CMGRs.csv')

genes_to_remove <- c("PLPP1", "PLPP2", "PLPP3")
filtered_gene <- gene[!gene$SYMBOL %in% genes_to_remove, ]
gene_set <- data.frame(symbols = filtered_gene, pathway =rep('score',length(filtered_gene)))
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
load("../00_rawdata/tcga_mrna.rdata")
library(GSVA)
expr <- tcga_mrna[, -(1:2)] 
rownames(expr) <- tcga_mrna$sample_id
expr <- t(expr)
ssgsea_res <-gsva(as.matrix(expr),
                  gene_list,
                  method = "ssgsea",
                  ssgsea.norm = TRUE, 
                  verbose = TRUE)

save(ssgsea_res, file = "ssgsea_score.rdata")
load("ssgsea_score.rdata")
ssgsea_df <- data.frame(ssgsea_res)%>%t(.)
ssgsea_df <-as.data.frame(ssgsea_df)
ssgsea_df$sample_id <- rownames(ssgsea_df)
ssgsea_df$sample_id <- gsub("\\.", "-", ssgsea_df$sample_id)

cancer_types <- tcga_mrna[, c("sample_id", "project")]
cancer_types <- unique(cancer_types)  


ssgsea_df <- merge(ssgsea_df, cancer_types, by = "sample_id")


ssgsea_long <- ssgsea_df %>%
  gather(key = "pathway", value = "ssgsea_score", -sample_id, -project)











p <- ggplot(ssgsea_long, aes(x = project, y = ssgsea_score, color = project)) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.7) +
  labs(x = "Cancer Types", y = "ssGSEA Score", title = " ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")


ggsave("01.ssgsea_scores_cancer_types.pdf", plot = p, width = 10, height = 6)
ggsave("01.ssgsea_scores_cancer_types.png", plot = p, width = 10, height = 6, dpi = 600)



library(GSVA)
library(GSEABase)
library(tidyverse)
gene_set1 <- getGmt("h.all.v2023.2.Hs.symbols.gmt")
ssgsea_res1 <-gsva(as.matrix(expr),
                   gene_set1,
                  method = "ssgsea",
                  ssgsea.norm = TRUE, 
                  verbose = TRUE)
save(ssgsea_res1, file = "ssgsea_score_Hallmark.rdata")

load("ssgsea_score_Hallmark.rdata")
ssgsea_df1 <- data.frame(ssgsea_res1) %>% t(.) %>% as.data.frame()
ssgsea_df1$sample_id <- rownames(ssgsea_df1)
ssgsea_df1$sample_id <- gsub("\\.", "-", ssgsea_df1$sample_id)


combined_df <- merge(ssgsea_df, ssgsea_df1, by = "sample_id")
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)


combined_long <- combined_df %>%
  gather(key = "pathway", value = "hallmark_score", -sample_id, -project, -score) %>%
  group_by(project, pathway) %>%
  summarise(correlation = cor(score, hallmark_score),
            p_value = cor.test(score, hallmark_score)$p.value)


combined_long <- combined_long %>%
  mutate(log_p_value = -log10(p_value),
         color = correlation)
write.csv(combined_long, file = '02.cor.csv', row.names = T, quote = F)



bubble_plot <- ggplot(combined_long, aes(x = project, y = pathway, size = log_p_value, color = correlation)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 5)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Cancer Type", y = "Pathway", color = "Correlation", size = "-log10(p-value)",
       title = " ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("02.correlation_bubble_plot.pdf", plot = bubble_plot, width = 12, height = 8)
ggsave("02.correlation_bubble_plot.png", plot = bubble_plot, width = 12, height = 8, dpi = 600)
