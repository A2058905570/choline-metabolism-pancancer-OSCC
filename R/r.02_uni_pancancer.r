
rm(list = ls())
setwd("./02_uncox_pan")


library(survival)
library(survminer)


OS <- read.csv('../00_rawdata/survival.tcga.csv')
head(OS)


OS1 <- OS[, c('sample_id', 'project', 'OS', 'OS.time')]
head(OS1)


expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
head(expr[, 1:6])


merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)
head(merged_data[, 1:6])


tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]


cox_results <- list()


for (proj in unique(merged_data$project)) {
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {
    

    cox_ph_test <- cox.zph(coxph(Surv(OS.time, OS) ~ get(gene), data = cox.dat))
    

    if (cox_ph_test[["table"]][1, "p"] >= 0.05) {
      cox_model <- coxph(Surv(OS.time, OS) ~ get(gene), data = cox.dat)
      cox_summary <- summary(cox_model)
      

      p.value <- signif(cox_summary$wald["pvalue"], digits = 3)
      HR <- signif(exp(cox_summary$coef[1]), digits = 4)  
      HR.confint.lower <- signif(exp(cox_summary$conf.int[,"lower .95"]), 4)
      HR.confint.upper <- signif(exp(cox_summary$conf.int[,"upper .95"]), 4)
      HR_with_CI <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
      

      result <- c(p.value, HR, HR_with_CI)
      names(result) <- c("p.value", "HR", "HR (95% CI for HR)")
      cox_results[[paste(proj, gene, sep = "_")]] <- c(proj, result)  
    } else {
      cat("Gene", gene, "does not satisfy PH assumption (p =", cox_ph_test[["table"]][1, "p"], "), skipping Cox regression.\n")
    }
  }
}


result_df <- data.frame(
  Project_Gene = character(),
  Project = character(),
  Gene = character(),
  P_Value = numeric(),
  HR = numeric(),
  HR_95_CI = character(),
  stringsAsFactors = FALSE
)


for (key in names(cox_results)) {

  project_gene <- key
  project <- strsplit(key, "_")[[1]][1]
  gene <- strsplit(key, "_")[[1]][2]
  

  p_value <- as.numeric(cox_results[[key]][["p.value"]])
  hr <- as.numeric(cox_results[[key]][["HR"]])
  hr_with_ci <- cox_results[[key]][["HR (95% CI for HR)"]]
  

  result_df <- rbind(result_df, data.frame(
    Project_Gene = project_gene,
    Project = project,
    Gene = gene,
    P_Value = p_value,
    HR = hr,
    HR_95_CI = hr_with_ci,
    stringsAsFactors = FALSE
  ))
}


str(result_df)



head(result_df)

result_file <- "uni_cox_regression_results.csv"
write.csv(result_df, file = result_file, row.names = FALSE)

colnames(result_df)

significant_results_df <- result_df[result_df$P_Value < 0.05, ]
write.csv(significant_results_df, "uin_sig_cox_results.csv", row.names = FALSE)



rm(list = ls())



setwd("./02_uncox_pan")


result_df <- read.csv('uni_cox_regression_results.csv')
head(result_df)
library(dplyr)
library(tidyr)
library(ggplot2)


result_df <- result_df %>%
  mutate(
    P_Value = as.numeric(P_Value),
    HR = as.numeric(HR),
    Color = case_when(
      P_Value < 0.05 & HR > 1 ~ "red",
      P_Value < 0.05 & HR < 1 ~ "blue",
      TRUE ~ "grey"
    )
  )



library(ComplexHeatmap)
library(circlize)  


heatmap_matrix <- result_df %>%
  dplyr::select(Project, Gene, Color) %>%
  tidyr::spread(key = Project, value = Color, fill = "grey")


color_mapping <- c("grey" = 0, "blue" = -1, "red" = 1)
heatmap_matrix_num <- heatmap_matrix %>%
  mutate(across(-Gene, ~ color_mapping[.])) %>%
  column_to_rownames(var = "Gene")


col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "grey", "red"))



ht<-Heatmap(as.matrix(heatmap_matrix_num),
        name = "Significance",
        col = col_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "Cancer Type",
        row_title = "Gene",
        cluster_columns = FALSE,  
        cluster_rows =   FALSE,
        heatmap_legend_param = list(
          title = "Significance",
          at = c(-1, 0, 1),
          labels = c("Protective", " ", "Risky")
        ))


pdf("01.heatmap.pdf", width = 8, height = 12)
draw(ht)
dev.off()


png("01.heatmap.png", width = 800, height = 1200, res = 100)
draw(ht)
dev.off()



rm(list = ls())


setwd("./02_uncox_pan")


result_df <- read.csv('uni_cox_regression_results.csv')

significant_results <- result_df[result_df$P_Value < 0.05, ]



gene_summary <- significant_results %>%
  mutate(HR_Category = ifelse(HR > 1, "HR > 1", "HR < 1")) %>%
  group_by(Gene, HR_Category) %>%
  summarise(Significant_Count = n()) %>%
  spread(HR_Category, Significant_Count, fill = 0) %>%
  mutate(Total = `HR > 1` + `HR < 1`)








gene_long <- gather(gene_summary, key = "HR_Category", value = "Significant_Count", -Gene, -Total)


















gene_long <- gene_long %>%
  mutate(Significant_Count = ifelse(HR_Category == "HR < 1", -Significant_Count, Significant_Count))


library(ggplot2)
library(tidyr)

p<-ggplot(gene_long, aes(x = reorder(Gene, Total), y = Significant_Count, fill = HR_Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  labs(
    title = "Number of Significant Genes (p < 0.05) Across All Cancer Types",
    x = "Gene",
    y = "Number of Cancer Types",
    fill = "HR Category"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("HR < 1" = "blue", "HR > 1" = "red")) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  

pdf("02.significant_genes_barplot.pdf", width = 8, height = 12)
print(p)
dev.off()


ggsave("02.significant_genes_barplot.png", plot = p, width = 8, height = 12, dpi = 600)
print(p)
dev.off()

















rm(list = ls())

setwd('./02_uncox_pan')
library(lance)
library(readxl)
library(readr)
library(tidyverse)
expr<-read.csv('../00_rawdata/tcga.tpm.csv',check.names = F,row.names = 1)
head(expr[,1:6])
id<-read.csv("../00_rawdata/OSCC.id.csv")
id<-id$x
expr1 <- expr[expr$sample_id%in%id, ]

col_names<-expr1$sample_id
tumor_sample_cols <- col_names[substr(col_names, 14, 14) == "0"]
expr1_subset <- expr1[expr1$sample_id %in% tumor_sample_cols, ]
head(expr1_subset[,1:6])

df <- expr1_subset[, -which(names(expr1_subset) == "project")]
head(df[,1:6])
write.csv(df, file = "OSCC_TPM.csv", quote = TRUE, row.names = FALSE)



OS<-read.csv('../00_rawdata/survival.tcga.csv')
head(OS)
OS1<-OS[,c('sample_id','OS','OS.time')]
head(OS1)

OS1_BRCA_subset <- OS1[OS1$sample_id %in% tumor_sample_cols, ]
write.csv(OS1_BRCA_subset, file = "OSCC_OS.csv", quote = TRUE, row.names = FALSE)

train_dat <- merge(OS1_BRCA_subset, df, by = 'sample_id')
train_dat <- column_to_rownames(train_dat, var = 'sample_id')
train_data <- train_dat %>% na.omit()


head(train_data[,1:6])


gene0<-colnames(df)
gene <- gene0[-1]
gene <- gsub('-','_',gene)
colnames(train_data) <- gsub('-','_',colnames(train_data))

library(survival)
res.all <- data.frame()
for (i in 1:length(gene)) {

  gene_list <- gene[i]
  cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                                gene_list))
  
  cox_more <- coxph(cox_data,
                    data = train_data)
  cox_zph <- cox.zph(cox_more)
  cox_table <- cox_zph$table[-nrow(cox_zph$table),]
  cox_table
  res <- data.frame(gene=gene[i],PH.p=cox_table[3])
  res.all <- rbind(res.all, res)
}

ph.gene <- res.all[as.numeric(res.all$PH.p)>0.05,]
ph.gene$gene <- gsub('_','-',ph.gene$gene)
colnames(train_data) <- gsub('_','-',colnames(train_data))
dim(ph.gene)

write.csv(res.all,
          file = "01.OSCC_uniCox_PH.csv",
          quote = F,
          row.names = F)


train_data <- cbind(train_data[,1:2],train_data[,ph.gene$gene])
library(survival)
library(survminer)
outTab <- data.frame() 
for(i in colnames(train_data[,3:ncol(train_data)])){
  cox <- coxph(Surv(OS.time,OS) ~ train_data[,i], data = train_data)
  coxSummary <- summary(cox)
  outTab <- rbind(outTab,
                  cbind(id=i,
                        HR=coxSummary$conf.int[,"exp(coef)"],
                        HR.95L=coxSummary$conf.int[,"lower .95"],
                        HR.95H=coxSummary$conf.int[,"upper .95"],
                        pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
head(outTab)
outTab <- outTab[order(as.numeric(outTab$HR), decreasing = T),]
write.csv(outTab, file="01.OSCC_uniCox_all.csv", row.names=F, quote=F)

res_results_0.05 <- outTab[which(as.numeric(outTab$pvalue) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05) 




res_results_0.05 <- res_results_0.05[order(as.numeric(res_results_0.05$HR.95H), decreasing = T),]



dim(res_results_0.05)


write.csv(res_results_0.05,
          file = "01.OSCC_uniCox_0.05.csv",
          quote = F,
          row.names = F)



rm(list = ls())



setwd("./02_uncox_pan")


library(survival)
library(survminer)


OS <- read.csv('../00_rawdata/survival.tcga.csv')
head(OS)


OS1 <- OS[, c('sample_id', 'project', 'DSS', 'DSS.time')]
head(OS1)


expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
head(expr[, 1:6])


merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)
head(merged_data[, 1:6])
merged_data <- merged_data[!is.na(merged_data$DSS.time), ]

tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]


cox_results <- list()


for (proj in unique(merged_data$project)) {
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {
    

    cox_ph_test <- cox.zph(coxph(Surv(DSS.time, DSS) ~ get(gene), data = cox.dat))
    

    if (cox_ph_test[["table"]][1, "p"] >= 0.05) {
      cox_model <- coxph(Surv(DSS.time, DSS) ~ get(gene), data = cox.dat)
      cox_summary <- summary(cox_model)
      

      p.value <- signif(cox_summary$wald["pvalue"], digits = 3)
      HR <- signif(exp(cox_summary$coef[1]), digits = 4)  
      HR.confint.lower <- signif(exp(cox_summary$conf.int[,"lower .95"]), 4)
      HR.confint.upper <- signif(exp(cox_summary$conf.int[,"upper .95"]), 4)
      HR_with_CI <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
      

      result <- c(p.value, HR, HR_with_CI)
      names(result) <- c("p.value", "HR", "HR (95% CI for HR)")
      cox_results[[paste(proj, gene, sep = "_")]] <- c(proj, result)  
    } else {
      cat("Gene", gene, "does not satisfy PH assumption (p =", cox_ph_test[["table"]][1, "p"], "), skipping Cox regression.\n")
    }
  }
}


result_df <- data.frame(
  Project_Gene = character(),
  Project = character(),
  Gene = character(),
  P_Value = numeric(),
  HR = numeric(),
  HR_95_CI = character(),
  stringsAsFactors = FALSE
)


for (key in names(cox_results)) {

  project_gene <- key
  project <- strsplit(key, "_")[[1]][1]
  gene <- strsplit(key, "_")[[1]][2]
  

  p_value <- as.numeric(cox_results[[key]][["p.value"]])
  hr <- as.numeric(cox_results[[key]][["HR"]])
  hr_with_ci <- cox_results[[key]][["HR (95% CI for HR)"]]
  

  result_df <- rbind(result_df, data.frame(
    Project_Gene = project_gene,
    Project = project,
    Gene = gene,
    P_Value = p_value,
    HR = hr,
    HR_95_CI = hr_with_ci,
    stringsAsFactors = FALSE
  ))
}


str(result_df)

head(result_df)

result_file <- "DSS_uni_cox_regression_results.csv"
write.csv(result_df, file = result_file, row.names = FALSE)




rm(list = ls())


setwd("./02_uncox_pan")


library(survival)
library(survminer)


OS <- read.csv('../00_rawdata/survival.tcga.csv')
head(OS)


OS1 <- OS[, c('sample_id', 'project', 'DFI', 'DFI.time')]
head(OS1)


expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
head(expr[, 1:6])


merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)
head(merged_data[, 1:6])
merged_data <- merged_data[!is.na(merged_data$DFI.time), ]

tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]


cox_results <- list()


for (proj in unique(merged_data$project)) {
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {
    

    cox_ph_test <- cox.zph(coxph(Surv(DFI.time, DFI) ~ get(gene), data = cox.dat))
    

    if (cox_ph_test[["table"]][1, "p"] >= 0.05) {
      cox_model <- coxph(Surv(DFI.time, DFI) ~ get(gene), data = cox.dat)
      cox_summary <- summary(cox_model)
      

      p.value <- signif(cox_summary$wald["pvalue"], digits = 3)
      HR <- signif(exp(cox_summary$coef[1]), digits = 4)  
      HR.confint.lower <- signif(exp(cox_summary$conf.int[,"lower .95"]), 4)
      HR.confint.upper <- signif(exp(cox_summary$conf.int[,"upper .95"]), 4)
      HR_with_CI <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
      

      result <- c(p.value, HR, HR_with_CI)
      names(result) <- c("p.value", "HR", "HR (95% CI for HR)")
      cox_results[[paste(proj, gene, sep = "_")]] <- c(proj, result)  
    } else {
      cat("Gene", gene, "does not satisfy PH assumption (p =", cox_ph_test[["table"]][1, "p"], "), skipping Cox regression.\n")
    }
  }
}


result_df <- data.frame(
  Project_Gene = character(),
  Project = character(),
  Gene = character(),
  P_Value = numeric(),
  HR = numeric(),
  HR_95_CI = character(),
  stringsAsFactors = FALSE
)


for (key in names(cox_results)) {

  project_gene <- key
  project <- strsplit(key, "_")[[1]][1]
  gene <- strsplit(key, "_")[[1]][2]
  

  p_value <- as.numeric(cox_results[[key]][["p.value"]])
  hr <- as.numeric(cox_results[[key]][["HR"]])
  hr_with_ci <- cox_results[[key]][["HR (95% CI for HR)"]]
  

  result_df <- rbind(result_df, data.frame(
    Project_Gene = project_gene,
    Project = project,
    Gene = gene,
    P_Value = p_value,
    HR = hr,
    HR_95_CI = hr_with_ci,
    stringsAsFactors = FALSE
  ))
}


str(result_df)

head(result_df)

result_file <- "DFI_uni_cox_regression_results.csv"
write.csv(result_df, file = result_file, row.names = FALSE)



rm(list = ls())


setwd("./02_uncox_pan")


library(survival)
library(survminer)


OS <- read.csv('../00_rawdata/survival.tcga.csv')
head(OS)


OS1 <- OS[, c('sample_id', 'project', 'PFI', 'PFI.time')]
head(OS1)


expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
head(expr[, 1:6])


merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)
head(merged_data[, 1:6])
merged_data <- merged_data[!is.na(merged_data$PFI.time), ]

tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]


cox_results <- list()


for (proj in unique(merged_data$project)) {
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {
    

    cox_ph_test <- cox.zph(coxph(Surv(PFI.time, PFI) ~ get(gene), data = cox.dat))
    

    if (cox_ph_test[["table"]][1, "p"] >= 0.05) {
      cox_model <- coxph(Surv(PFI.time, PFI) ~ get(gene), data = cox.dat)
      cox_summary <- summary(cox_model)
      

      p.value <- signif(cox_summary$wald["pvalue"], digits = 3)
      HR <- signif(exp(cox_summary$coef[1]), digits = 4)  
      HR.confint.lower <- signif(exp(cox_summary$conf.int[,"lower .95"]), 4)
      HR.confint.upper <- signif(exp(cox_summary$conf.int[,"upper .95"]), 4)
      HR_with_CI <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
      

      result <- c(p.value, HR, HR_with_CI)
      names(result) <- c("p.value", "HR", "HR (95% CI for HR)")
      cox_results[[paste(proj, gene, sep = "_")]] <- c(proj, result)  
    } else {
      cat("Gene", gene, "does not satisfy PH assumption (p =", cox_ph_test[["table"]][1, "p"], "), skipping Cox regression.\n")
    }
  }
}


result_df <- data.frame(
  Project_Gene = character(),
  Project = character(),
  Gene = character(),
  P_Value = numeric(),
  HR = numeric(),
  HR_95_CI = character(),
  stringsAsFactors = FALSE
)


for (key in names(cox_results)) {

  project_gene <- key
  project <- strsplit(key, "_")[[1]][1]
  gene <- strsplit(key, "_")[[1]][2]
  

  p_value <- as.numeric(cox_results[[key]][["p.value"]])
  hr <- as.numeric(cox_results[[key]][["HR"]])
  hr_with_ci <- cox_results[[key]][["HR (95% CI for HR)"]]
  

  result_df <- rbind(result_df, data.frame(
    Project_Gene = project_gene,
    Project = project,
    Gene = gene,
    P_Value = p_value,
    HR = hr,
    HR_95_CI = hr_with_ci,
    stringsAsFactors = FALSE
  ))
}


str(result_df)

head(result_df)

result_file <- "PFI_uni_cox_regression_results.csv"
write.csv(result_df, file = result_file, row.names = FALSE)

