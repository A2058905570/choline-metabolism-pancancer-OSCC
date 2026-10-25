rm(list = ls())

setwd('./04_OSCC_uniCox')
library(lance)
library(readxl)
library(readr)
library(tidyverse)
expr<-read.csv('../00_rawdata/tcga.tpm.csv',check.names = F,row.names = 1)
head(expr[,1:6])
id<-read.csv('../00_rawdata/OSCC.id.csv')
id<-id$x
expr1 <- expr[expr$sample_id%in% id, ]

col_names<-expr1$sample_id
tumor_sample_cols <- col_names[nchar(col_names) >= 14 & substr(col_names, 14, 14) == "0"]
expr1_subset <- expr1[expr1$sample_id %in% tumor_sample_cols, ]
head(expr1_subset[,1:6])

df <- expr1_subset[, -which(names(expr1_subset) == "project")]
head(df[,1:6])


OS<-read.csv('../00_rawdata/survival.tcga.csv')
head(OS)
OS1<-OS[,c('sample_id','OS','OS.time')]
head(OS1)
OS1_BRCA_subset <- OS1[OS1$sample_id %in% tumor_sample_cols, ]
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




res_results_plot <- res_results_0.05
res_results_plot <- res_results_plot[order(as.numeric(res_results_plot$pvalue), decreasing = F),]
res_results_plot <- head(res_results_plot,100)
res_results_plot <- res_results_plot[order(as.numeric(res_results_plot$HR), decreasing = F),]
res_results_plot$HR <- as.numeric(res_results_plot$HR)
res_results_plot$HR.95L <- as.numeric(res_results_plot$HR.95L)
res_results_plot$HR.95H <- as.numeric(res_results_plot$HR.95H)
res_results_plot$pvalue <- as.numeric(res_results_plot$pvalue)

hz <- paste(round(res_results_plot$HR,3),
            "(",round(res_results_plot$HR.95L,3),
            "-",round(res_results_plot$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA, "Gene",res_results_plot$id),
                   c(NA,"P value",ifelse(res_results_plot$pvalue<0.001,
                                         "< 0.001",
                                         round(res_results_plot$pvalue,4))),
                   c(NA,"Hazard Ratio(95% CI)",hz))








nrow(tabletext) + 1

library(forestplot)
pdf(file = "01.os_univariate_cox_forest.pdf", height = 15, width = 11, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(FALSE,TRUE,rep(FALSE, nrow(tabletext) + 1)),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,NA,res_results_plot$HR),
           lower=c(NA,NA,res_results_plot$HR.95L), 
           upper=c(NA,NA,res_results_plot$HR.95H), 
           boxsize=0.15,lwd.ci=2,   
           ci.vertices.height = 0.05,ci.vertices=TRUE, 
           zero=1,lwd.zero=2,    
           colgap=unit(5,"mm"),    
           xticks = c(0,1,5,9), 
           lwd.xaxis=2,            
           lineheight = unit(0.7,"cm"), 
           graphwidth = unit(.5,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.4, fontfamily = "Times"),
                          ticks=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.3, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.5, fontface = "bold", fontfamily="Times")),
           xlab="OS Hazard Ratio",
           grid = T, 
           clip = c(0,15)
)  
dev.off()


png(file = "01.os_univariate_cox_forest.png", height = 15, width = 11, units='in',res=600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(FALSE,TRUE,rep(FALSE, nrow(tabletext) + 1)),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,NA,res_results_plot$HR),
           lower=c(NA,NA,res_results_plot$HR.95L), 
           upper=c(NA,NA,res_results_plot$HR.95H), 
           boxsize=0.15,lwd.ci=2,   
           ci.vertices.height = 0.05,ci.vertices=TRUE, 
           zero=1,lwd.zero=2,    
           colgap=unit(5,"mm"),    
           xticks = c(0,1,5,9), 
           lwd.xaxis=2,            
           lineheight = unit(0.7,"cm"), 
           graphwidth = unit(.5,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.4, fontfamily = "Times"),
                          ticks=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.3, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.5, fontface = "bold", fontfamily="Times")),
           xlab="OS Hazard Ratio",
           grid = T, 
           clip = c(0,15)
)  
dev.off()



DEG1 <- read.csv("../01_CMGRs_expression/diff_OSCC.csv", row.names = 1)
DEG2 <- read.csv("01.OSCC_uniCox_0.05.csv", row.names = 1)
DEG3 <- read.csv("../03_km_pan/OSCC_sig_km_results.csv")
DEERG <- data.frame(symbol=intersect(rownames(DEG1),intersect(rownames(DEG2),DEG3$gene)))
write.csv(DEERG,file = '01.OSCC_gene.csv',row.names = F)


library(ggvenn)
mydata<-list('DEGs'=rownames(DEG1),'Uni_cox'=rownames(DEG2),'K-M'=DEG3$gene)
pdf('02.OSCC_venn.pdf',w=7,h=7,family='Times')
ggvenn(mydata,c('DEGs', 'Uni_cox', 'K-M'),
       fill_color = c("#49afca", "#d64c39","#fed53b"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 5,
       stroke_color="white",
       stroke_linetype="dashed",
       set_name_color='black',
       text_color = 'black')
dev.off()
png('02.OSCC_venn.png',w=7,h=7, units = 'in',res = 600,family='Times')
ggvenn(mydata,c('DEGs', 'Uni_cox', 'K-M'),
       fill_color = c("#49afca", "#d64c39","#fed53b"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 1,
       text_size = 5,
       stroke_color="white",
       stroke_linetype="dashed",
       set_name_color='black',
       text_color = 'black')
dev.off()





library(survival)
library(survminer)
cox_gene <- rownames(DEG2)
a<-c("OS","OS.time",cox_gene)
train_dat_cox_gene_sv<-train_data[,colnames(train_data)%in%a]


res.all.cox.gene <- data.frame()

for (i in 1:length(cox_gene)) {

  gene_list <- cox_gene[i]
  cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                                paste(gene_list,
                                      sep = '',
                                      collapse = '+')))
  
  cox_more <- coxph(cox_data,
                    data = train_dat_cox_gene_sv)
  cox_zph <- cox.zph(cox_more)
  cox_table <- cox_zph$table[-nrow(cox_zph$table),]
  cox_table
  res <- data.frame(gene=cox_gene[i],PH.p=cox_table[3])
  res.all <- rbind(res.all, res)
  p <- ggcoxzph(cox_zph)
  pdf(paste0('03.ph.test_',gene_list,'.pdf'), width = 6,height = 4,onefile = F,family='Times')
  print(p)
  dev.off()
  
  png(paste0('03.ph.test_',gene_list,'.png'), width = 6,height = 4,family='Times',units='in',res=600)
  print(p)
  dev.off()
}


