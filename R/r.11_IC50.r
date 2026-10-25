rm(list = ls())

setwd("./11_IC50")
tcga_mrna<-read.csv("../00_rawdata/output_mRNA_lncRNA_expr/TCGA-HNSC_mrna_expr_tpm.csv",check.names = F,row.names = 1)
new_colnames <- substr(colnames(tcga_mrna), 1, 15) 
colnames(tcga_mrna) <- new_colnames  
head(tcga_mrna[,1:4],4)
riskScore <- read.csv("../08_GSEA/risk.csv")
dat <- tcga_mrna[,riskScore$id]











library(pRRophetic)
library(ggplot2)
set.seed(12345)
drug<-read.table(file = '/data/nas3/mayanyu/data/drugs.txt',sep='\t',header=F)
a <- data.frame(row.names=riskScore$id,riskScore$id)
cnt<-1
while (cnt < 139) {
  predictedPtype <- pRRopheticPredict(as.matrix(dat), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
write.csv(a,'IC50.csv')

a <- read.csv("IC50.csv",header = TRUE,row.names = 1)
colnames(a)[1] <- "id"

result <- a[,-1]

removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
result<-removeColsAllNa(result)
na_flag <- apply(is.na(result), 2, sum)
result <- result[, which(na_flag == 0)]
dim(result)

medicinal_result <- t(result)

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

high_group <- riskScore$id[which(riskScore$risk =="High Risk")]
low_group <- riskScore$id[which(riskScore$risk =="Low Risk")]


pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (j in 1:nrow(medicinal_result)){
  pvalue[j, 1] = p.value = wilcox.test(medicinal_result[j, high_group],
                                       medicinal_result[j, low_group])$p.value
  log2FoldChange[j, 1] = mean(medicinal_result[j, high_group]) - 
    mean(medicinal_result[j, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****", sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)
write.csv(rTable, file = "02.ic50_wilcox_test.csv")
DE.drug<-rTable[which(rTable$pvalue<0.05),]
DE.drug_up <- DE.drug[which(DE.drug$log2FoldChange>0),]
DE.drug_down <- DE.drug[which(DE.drug$log2FoldChange<0),]
write.csv(DE.drug_up, file = "03.ic50_up.csv")
write.csv(DE.drug_down, file = "04.ic50_down.csv")

all(rownames(rTable) == rownames(medicinal_result))

DE.drug_down <- DE.drug_down[order(DE.drug_down$padj),]
DE.drug_down <- DE.drug_down[c(1:5),]
if (length(rownames(DE.drug_down)) != 0) {
  drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
  drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_down),] 
  drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
  violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
  violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                              "High Risk", "Low Risk") 
  violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High Risk", "Low Risk"))

  drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                  color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                  add = "jitter",
                                  short.panel.labs = T,
                                  ggtheme = theme_bw()) +
    stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
  drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                             facet.by = "drugs",
                             short.panel.labs = T,
                             panel.labs.background = list(fill = "white"),
                             ncol = 5,
                             scales = "free_y") + xlab("") + ylab("IC(50)") +


    theme(panel.grid = element_blank(),
          legend.position = "none",

          axis.title.y = element_text(size = 18, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          text = element_text(size = 13, face = "bold"))
  drugs_hub_boxplot
  ggsave(filename = "ic50_drugs.plot(down).pdf", height = 5, width = 12)
  ggsave(filename = "ic50_drugs.plot(down).png", height = 5, width = 12)
}



DE.drug_up <- DE.drug_up[order(DE.drug_up$padj),]
DE.drug_up <- DE.drug_up[c(1:5),]
if (length(rownames(DE.drug_up)) != 0){
  drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
  drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_up),]

  drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
  violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
  violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                              "High Risk", "Low Risk") 
  violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High Risk", "Low Risk"))
  head(violin_dat)
  drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                  color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                  add = "jitter",
                                  short.panel.labs = T,
                                  ggtheme = theme_bw()) +
    stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
  drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                             facet.by = "drugs",
                             short.panel.labs = T,
                             panel.labs.background = list(fill = "white"),
                             ncol = 5,
                             scales = "free_y") + xlab("") + ylab("IC(50)") +


    theme(panel.grid = element_blank(),
          legend.position = "none",

          axis.title.y = element_text(size = 18, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          text = element_text(size = 13, face = "bold"))
  drugs_hub_boxplot
  ggsave(filename = "ic50_drugs.plot(up).pdf", height = 5, width = 12)
  ggsave(filename = "ic50_drugs.plot(up).png", height = 5, width = 12)
}








































































dat_up<-result
dat_score<- read.csv('../05_model/01.OSCC_train_risk.csv')
riskscore<-dat_score[,"riskScore"]
riskscore<-as.data.frame(riskscore)
rownames(riskscore)<-dat_score$X


sorted_indices <- order(riskscore$riskscore)
dat_up <- dat_up[sorted_indices, ]
riskscore <- riskscore[sorted_indices, , drop = FALSE]


cor_results <- apply(dat_up, 2, function(col) {
  cor.test(col, riskscore$riskscore, method = "spearman")
})

cor_values <- sapply(cor_results, function(x) x$estimate)
p_values <- sapply(cor_results, function(x) x$p.value)
Pathway_names <- names(dat_up)  


cor_data <- data.frame(
  Pathway = Pathway_names,
  SpearmanR = cor_values,
  PValue = p_values
)
write.csv(cor_data,file = "all_cor_data.csv")




sorted_cor_data <- cor_data[order(abs(cor_data$SpearmanR), decreasing = TRUE), ]  


pos_top5 <- sorted_cor_data[sorted_cor_data$SpearmanR > 0, ][1:5, ]  


neg_top5 <- sorted_cor_data[sorted_cor_data$SpearmanR < 0, ][1:5, ]  


plot_data <- rbind(pos_top5, neg_top5)  
library(ggplot2)  


ggplot(plot_data, aes(x = SpearmanR, y = reorder(Pathway, SpearmanR))) +  
  geom_segment(aes(x = 0, xend = SpearmanR, y = reorder(Pathway, SpearmanR), yend = reorder(Pathway, SpearmanR)), color = "grey") +  
  geom_point(aes(size = -log10(PValue), color = SpearmanR), show.legend = TRUE) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
  labs(x = "Spearman Correlation", y = " ", title = " ") +  
  theme_minimal() +  
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  
  scale_size_continuous(name = "-log10(P-value)", range = c(2, 10))  
ggsave("spearman_correlation.pdf", width = 6, height = 8)
ggsave("spearman_correlation.png", width = 6, height = 8)

