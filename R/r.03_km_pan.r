
rm(list = ls())

setwd("./03_km_pan")


OS <- read.csv('../00_rawdata/survival.tcga.csv')

expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')

OS1 <- OS[, c('sample_id', 'project', 'OS', 'OS.time')]
merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)

tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]
library(survival)   

cancer_types <- unique(merged_data$cancer_type)  

library(survival)
library(survminer)

km_results <- list()


km_results_df <- data.frame(project = character(),
                            gene = character(),
                            cutpoint = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)


for (proj in unique(merged_data$project)) {  
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]  
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {  

    tryCatch({
      cutpoint <- surv_cutpoint(cox.dat, time = "OS.time", event = "OS", variables = gene)
      best_cut <- cutpoint$cutpoint[[1]]

      cox.dat$expr_group <- ifelse(cox.dat[[gene]] > best_cut, "High", "Low")

      surv_obj <- Surv(time = cox.dat$OS.time, event = cox.dat$OS)
      km_fit <- survfit(surv_obj ~ expr_group, data = cox.dat)

      surv_diff <- survdiff(surv_obj ~ expr_group, data = cox.dat)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      

      km_results_df <- rbind(km_results_df, data.frame(
        project = proj,
        gene = gene,
        cutpoint = best_cut,
        pvalue = p_value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat("Error in gene:", gene, "of project:", proj, "\n")
      cat("Error message:", e$message, "\n")
    })
  }  
}


write.csv(km_results_df, "cut_all_km_results.csv", row.names = T)

significant_results_df <- km_results_df[km_results_df$pvalue < 0.05, ]

write.csv(significant_results_df, "cut_sig_km_results.csv", row.names = T)



rm(list = ls())

setwd("./03_km_pan")

km_results_df <- read.csv('cut_all_km_results.csv')
HR<-read.csv('../02_uncox_pan/uni_cox_regression_results.csv')
colnames(HR)

km_results_df <- km_results_df %>%
  mutate(Project_Gene = paste(project, gene, sep = "_"))
merged_df <- inner_join(km_results_df, HR, by = "Project_Gene")

km_results_df<-merged_df[,c("project","gene","cutpoint","pvalue","HR")]




km_results_df$neg_log10_pvalue <- -log10(km_results_df$pvalue)

km_results_df <- km_results_df %>%
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "Not Significant"))
p <- ggplot(km_results_df, aes(x = project, y = gene, size = neg_log10_pvalue, color = HR, shape = significance)) +
  geom_point(alpha = 1) +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 5)) +
  scale_color_gradient2(name = "HR", low = "darkblue", mid = "white", high = "#d73027", midpoint = 1) +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Not Significant" = 1)) +
  labs(
    x = "Project",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
p

ggsave("01.OS.pdf", plot = p, width = 8, height = 8)


ggsave("01.OS.png", plot = p, width = 8, height = 8, dpi = 600)




  






OS <- read.csv('../00_rawdata/survival.OSCC.csv')

expr <- read.csv('../00_rawdata/OSCC.tpm.csv')
expr<-expr[, -which(names(expr) == "project")]

OS1 <- OS[, c('sample_id', 'project', 'OS', 'OS.time')]
merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)

tumor_samples <- merged_data$sample_id[substr(merged_data$sample_id, 14, 14) == "0"]
library(survival)   

cancer_types <- unique(merged_data$project)  

library(survival)
library(survminer)

km_results <- list()


km_results_df <- data.frame(project = character(),
                            gene = character(),
                            cutpoint = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)


for (proj in unique(merged_data$project)) {  
  cox.dat <- merged_data[merged_data$project == proj & merged_data$sample_id %in% tumor_samples, ]  
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {  

    tryCatch({
      cutpoint <- surv_cutpoint(cox.dat, time = "OS.time", event = "OS", variables = gene)
      best_cut <- cutpoint$cutpoint[[1]]

      cox.dat$expr_group <- ifelse(cox.dat[[gene]] > best_cut, "High", "Low")

      surv_obj <- Surv(time = cox.dat$OS.time, event = cox.dat$OS)
      km_fit <- survfit(surv_obj ~ expr_group, data = cox.dat)

      surv_diff <- survdiff(surv_obj ~ expr_group, data = cox.dat)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      

      km_results_df <- rbind(km_results_df, data.frame(
        project = proj,
        gene = gene,
        cutpoint = best_cut,
        pvalue = p_value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat("Error in gene:", gene, "of project:", proj, "\n")
      cat("Error message:", e$message, "\n")
    })
  }  
}


write.csv(km_results_df, "OSCC_cut_all_km_results.csv", row.names = FALSE)

significant_results_df <- km_results_df[km_results_df$pvalue < 0.05, ]

write.csv(significant_results_df, "OSCC_sig_km_results.csv", row.names = FALSE)






setwd("./03_km_pan")
OS <- read.csv('../00_rawdata/survival.tcga.csv')
expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
colnames(OS)
OS1 <- OS[, c('sample_id', 'project', 'DFI', 'DFI.time')]
merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)

merged_data$DFI.time <- as.numeric(as.character(merged_data$DFI.time))


filtered_data <- merged_data[!is.na(merged_data$DFI.time), ]


tumor_samples <- filtered_data$sample_id[substr(filtered_data$sample_id, 14, 14) == "0"]
library(survival)   

cancer_types <- unique(filtered_data$cancer_type)  

library(survival)
library(survminer)

km_results <- list()


km_results_df <- data.frame(project = character(),
                            gene = character(),
                            cutpoint = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)


for (proj in unique(filtered_data$project)) {  
  cox.dat <- filtered_data[filtered_data$project == proj & filtered_data$sample_id %in% tumor_samples, ]  
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {  

    tryCatch({
      cutpoint <- surv_cutpoint(cox.dat, time = "DFI.time", event = "DFI", variables = gene)
      best_cut <- cutpoint$cutpoint[[1]]

      cox.dat$expr_group <- ifelse(cox.dat[[gene]] > best_cut, "High", "Low")

      surv_obj <- Surv(time = cox.dat$DFI.time, event = cox.dat$DFI)
      km_fit <- survfit(surv_obj ~ expr_group, data = cox.dat)

      surv_diff <- survdiff(surv_obj ~ expr_group, data = cox.dat)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      

      km_results_df <- rbind(km_results_df, data.frame(
        project = proj,
        gene = gene,
        cutpoint = best_cut,
        pvalue = p_value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat("Error in gene:", gene, "of project:", proj, "\n")
      cat("Error message:", e$message, "\n")
    })
  }  
}


print(km_results_df)

write.csv(km_results_df, "DFI_cut_all_km_results.csv", row.names = T)

km_results_df <- read.csv('DFI_cut_all_km_results.csv')
HR<-read.csv('../02_uncox_pan/DFI_uni_cox_regression_results.csv')
colnames(HR)

km_results_df <- km_results_df %>%
  mutate(Project_Gene = paste(project, gene, sep = "_"))
merged_df <- inner_join(km_results_df, HR, by = "Project_Gene")

km_results_df<-merged_df[,c("project","gene","cutpoint","pvalue","HR")]




km_results_df$neg_log10_pvalue <- -log10(km_results_df$pvalue)

km_results_df <- km_results_df %>%
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "Not Significant"))
p <- ggplot(km_results_df, aes(x = project, y = gene, size = neg_log10_pvalue, color = HR, shape = significance)) +
  geom_point(alpha = 1) +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 5)) +
  scale_color_gradient2(name = "HR", low = "darkblue", mid = "white", high = "#d73027", midpoint = 1) +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Not Significant" = 1)) +
  labs(
    #title = "Kaplan-Meier Survival Analysis Results",
    x = "Project",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("02.DFI.pdf", plot = p, width = 8, height = 8)


ggsave("02.DFI.png", plot = p, width = 8, height = 8, dpi = 600)







rm(list = ls())

setwd("./03_km_pan")
OS <- read.csv('../00_rawdata/survival.tcga.csv')
expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
colnames(OS)
OS1 <- OS[, c('sample_id', 'project', 'DSS', 'DSS.time')]
merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)


merged_data$DSS.time <- as.numeric(as.character(merged_data$DSS.time))


filtered_data <- merged_data[!is.na(merged_data$DSS.time), ]


tumor_samples <- filtered_data$sample_id[substr(filtered_data$sample_id, 14, 14) == "0"]
library(survival)   

cancer_types <- unique(filtered_data$cancer_type)  

library(survival)
library(survminer)

km_results <- list()


km_results_df <- data.frame(project = character(),
                            gene = character(),
                            cutpoint = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)


for (proj in unique(filtered_data$project)) {  
  cox.dat <- filtered_data[filtered_data$project == proj & filtered_data$sample_id %in% tumor_samples, ]  
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {  

    tryCatch({
      cutpoint <- surv_cutpoint(cox.dat, time = "DSS.time", event = "DSS", variables = gene)
      best_cut <- cutpoint$cutpoint[[1]]

      cox.dat$expr_group <- ifelse(cox.dat[[gene]] > best_cut, "High", "Low")

      surv_obj <- Surv(time = cox.dat$DSS.time, event = cox.dat$DSS)
      km_fit <- survfit(surv_obj ~ expr_group, data = cox.dat)

      surv_diff <- survdiff(surv_obj ~ expr_group, data = cox.dat)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      

      km_results_df <- rbind(km_results_df, data.frame(
        project = proj,
        gene = gene,
        cutpoint = best_cut,
        pvalue = p_value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat("Error in gene:", gene, "of project:", proj, "\n")
      cat("Error message:", e$message, "\n")
    })
  }  
}


print(km_results_df)


write.csv(km_results_df, "DSS_cut_all_km_results.csv", row.names = T)

significant_results_df <- km_results_df[km_results_df$pvalue < 0.05, ]

write.csv(significant_results_df, "DSS_cut_sig_km_results.csv", row.names = T)

km_results_df <- read.csv('DSS_cut_all_km_results.csv')
HR<-read.csv('../02_uncox_pan/DSS_uni_cox_regression_results.csv')
colnames(HR)

km_results_df <- km_results_df %>%
  mutate(Project_Gene = paste(project, gene, sep = "_"))
merged_df <- inner_join(km_results_df, HR, by = "Project_Gene")

km_results_df<-merged_df[,c("project","gene","cutpoint","pvalue","HR")]




km_results_df$neg_log10_pvalue <- -log10(km_results_df$pvalue)

km_results_df <- km_results_df %>%
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "Not Significant"))
p <- ggplot(km_results_df, aes(x = project, y = gene, size = neg_log10_pvalue, color = HR, shape = significance)) +
  geom_point(alpha = 1) +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 5)) +
  scale_color_gradient2(name = "HR", low = "darkblue", mid = "white", high = "#d73027", midpoint = 1) +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Not Significant" = 1)) +
  labs(
    #title = "Kaplan-Meier Survival Analysis Results",
    x = "Project",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
p

ggsave("02.DSS.pdf", plot = p, width = 8, height = 8)


ggsave("02.DSS.png", plot = p, width = 8, height = 8, dpi = 600)





rm(list = ls())

setwd("./03_km_pan")
OS <- read.csv('../00_rawdata/survival.tcga.csv')
expr <- read.csv('../01_CMGRs_expression/dat_expr.csv')
colnames(OS)
OS1 <- OS[, c('sample_id', 'project', 'PFI', 'PFI.time')]
merged_data <- merge(OS1, expr, by = "sample_id", all.x = TRUE, sort = FALSE)
merged_data$PFI.time <- as.numeric(as.character(merged_data$PFI.time))
filtered_data <- merged_data[!is.na(merged_data$PFI.time), ]

tumor_samples <- filtered_data$sample_id[substr(filtered_data$sample_id, 14, 14) == "0"]
library(survival)   

cancer_types <- unique(filtered_data$cancer_type)  

library(survival)
library(survminer)

km_results <- list()


km_results_df <- data.frame(project = character(),
                            gene = character(),
                            cutpoint = numeric(),
                            pvalue = numeric(),
                            stringsAsFactors = FALSE)


for (proj in unique(filtered_data$project)) {  
  cox.dat <- filtered_data[filtered_data$project == proj & filtered_data$sample_id %in% tumor_samples, ]  
  

  genes <- colnames(cox.dat)[5:ncol(cox.dat)]  
  

  for (gene in genes) {  

    tryCatch({
      cutpoint <- surv_cutpoint(cox.dat, time = "PFI.time", event = "PFI", variables = gene)
      best_cut <- cutpoint$cutpoint[[1]]

      cox.dat$expr_group <- ifelse(cox.dat[[gene]] > best_cut, "High", "Low")

      surv_obj <- Surv(time = cox.dat$PFI.time, event = cox.dat$PFI)
      km_fit <- survfit(surv_obj ~ expr_group, data = cox.dat)

      surv_diff <- survdiff(surv_obj ~ expr_group, data = cox.dat)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      

      km_results_df <- rbind(km_results_df, data.frame(
        project = proj,
        gene = gene,
        cutpoint = best_cut,
        pvalue = p_value,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      cat("Error in gene:", gene, "of project:", proj, "\n")
      cat("Error message:", e$message, "\n")
    })
  }  
}



write.csv(km_results_df, "PFI_cut_all_km_results.csv", row.names = T)

km_results_df <- read.csv('PFI_cut_all_km_results.csv')
HR<-read.csv('../02_uncox_pan/PFI_uni_cox_regression_results.csv')
colnames(HR)

km_results_df <- km_results_df %>%
  mutate(Project_Gene = paste(project, gene, sep = "_"))
merged_df <- inner_join(km_results_df, HR, by = "Project_Gene")

km_results_df<-merged_df[,c("project","gene","cutpoint","pvalue","HR")]




km_results_df$neg_log10_pvalue <- -log10(km_results_df$pvalue)

km_results_df <- km_results_df %>%
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "Not Significant"))
p <- ggplot(km_results_df, aes(x = project, y = gene, size = neg_log10_pvalue, color = HR, shape = significance)) +
  geom_point(alpha = 1) +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 5)) +
  scale_color_gradient2(name = "HR", low = "darkblue", mid = "white", high = "#d73027", midpoint = 1) +
  scale_shape_manual(name = "Significance", values = c("Significant" = 16, "Not Significant" = 1)) +
  labs(
    #title = "Kaplan-Meier Survival Analysis Results",
    x = "Project",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
p

ggsave("02.PFI.pdf", plot = p, width = 8, height = 8)


ggsave("02.PFI.png", plot = p, width = 8, height = 8, dpi = 600)






