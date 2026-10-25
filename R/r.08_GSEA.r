rm(list = ls())

setwd("./08_GSEA")

library(DESeq2)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)
library(lance)
library(patchwork)

tcga_mrna<-read.csv("../00_rawdata/output_mRNA_lncRNA_expr/TCGA-HNSC_mrna_expr_count.csv",check.names = F,row.names = 1)
new_colnames <- substr(colnames(tcga_mrna), 1, 15) 
colnames(tcga_mrna) <- new_colnames  
head(tcga_mrna[,1:4],4)




dat_score<- read.csv('../05_model/01.OSCC_train_risk.csv') %>% lc.tableToNum()
dat_score$risk <- ifelse(dat_score$riskScore > 39.69651, 1, 0)
dat_score$risk <- ifelse(dat_score$risk == 1, 'High Risk', 'Low Risk')
colnames(dat_score)
risk <- dat_score[,c("X","risk")]
colnames(risk)<-c("id","risk")
write.csv(risk,file="risk.csv",row.names = F)
dat <- tcga_mrna[,risk$id]
dat_final <- dat
head(dat_final[,1:4],4)


group <- risk
colnames(group) <- c("sample","group")
group$group = factor(group$group, levels = c('High Risk','Low Risk'))
design.mat = cbind(Low = ifelse(group$group == "Low Risk", 1, 0), 
                   High = ifelse(group$group == "Low Risk", 0, 1))

contrast.mat = makeContrasts(contrasts="High-Low", levels=design.mat)
fit = lmFit(dat_final, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG=na.omit(fit)


risk_gene_symbol <- rownames(DEG)
risk_gene <- DEG  

gene <- bitr(risk_gene_symbol,
             fromType = "SYMBOL",
             toType = c("ENTREZID"),
             OrgDb = "org.Hs.eg.db")  
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)  
colnames(gene) <- c("GeneSymbol","ENTREZID")  
risk_gene$GeneSymbol <- rownames(risk_gene)  
data_all <- risk_gene %>% dplyr::inner_join(.,gene,by="GeneSymbol")    
data_all_sort <- data_all %>% arrange(desc(logFC))  
genelist <- data_all_sort$logFC 
names(genelist) <- data_all_sort$GeneSymbol 
write.csv(genelist,file="genelist.csv")

if(T){
  gsea.plot = function(res.kegg, top.hall,name){
    gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
    gsdata$Description = factor(gsdata$Description, levels = top.hall)
    p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + 

      theme_classic(14) +
      theme(panel.grid.major = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(colour = "grey92"), 
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_brewer(palette = "Paired") +
      ggtitle(paste0("Gene Set Enrichment Analysis: ",name)) +

      geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
      theme(legend.position = "top", 
            legend.justification = c(0,1),
            legend.title = element_blank(), 
            legend.background = element_rect(fill = "transparent")) +
      guides(
        color = guide_legend(
          ncol = 1, 
          byrow = TRUE,
          reverse = T)
      )+
      ylab("Running Enrichment Score") + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            text = element_text(face = "bold", family = "Times"))+
      theme(
        axis.title.y =element_text(size=15,family = "Times", face = "bold"),
        axis.text.y=element_text(size=15,family = "Times", face = "bold"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank()
      )
    i = 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1 }
    p2 = ggplot(gsdata, aes_(x = ~x)) + 
      geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + 
      xlab(NULL) + ylab(NULL) + 

      theme_classic(14) +
      theme(legend.position = "none", 
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank(),
            text = element_text(face = "bold", family = "Times")) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_color_brewer(palette = "Paired")+
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),

        legend.text=element_text(size=14))
    p3 <- gseaplot2(res.kegg,top.hall,subplots = 3)+
      
      theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
            axis.text.x =element_text(angle=0,size=15,hjust = 0,family = "Times", face = "bold"),
            axis.title.y =element_text(size=15,family = "Times", face = "bold"),
            axis.text.y=element_text(size=15,family = "Times", face = "bold"))+
      theme(
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.title=element_text(size=15) , legend.text=element_text(size=14))

    p <- p1/p2/p3+plot_layout(ncol = 1, height = c(4,1,2))
    return(p)
  }
}




kegggmt <- read.gmt("/data/nas3/mayanyu/data/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
KEGG<-GSEA(genelist,TERM2GENE = kegggmt,pvalueCutoff = 0.05) 
sortKEGG <- KEGG[KEGG$p.adjust<0.05 & abs(KEGG$NES)>1,]
dim(sortKEGG)
sortKEGG<-sortKEGG[order(sortKEGG$enrichmentScore, decreasing = F),]
write.csv(sortKEGG,file="KEGG_result.csv")
paths <- rbind(head(sortKEGG,5))

select_path <- paths$Description


pal <- gsea.plot(KEGG, select_path,"KEGG")
ggsave("KEGG.png",width = 11, height = 9,pal)
ggsave("KEGG.pdf",width = 11, height = 9,pal)















