
rm(list = ls())

setwd("./09_ssGSEA")

library(lance)

tcga_mrna<-read.csv("../00_rawdata/output_mRNA_lncRNA_expr/TCGA-HNSC_mrna_expr_tpm.csv",check.names = F,row.names = 1)
new_colnames <- substr(colnames(tcga_mrna), 1, 15) 
colnames(tcga_mrna) <- new_colnames  



group<-read.csv('../08_GSEA/risk.csv')
expr <- tcga_mrna[,group$id]
colnames(group) <- c('sample', 'group')
table(group$group)
case.sample<-group$sample[which(group$group=='High Risk')]
control.sample<-group$sample[which(group$group=='Low Risk')]
expr<-expr[,group$sample]

library(GSVA)
gene_set <- read.table("/data/nas3/mayanyu/data/mmc3.txt",
                       header = T,
                       sep ="\t")

dat.final <- as.matrix(expr)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.csv(ssgsea_score,
          file = "01.ssgsea_result_cell.csv",
          quote = F)


rm(list = ls())

setwd("./09_ssGSEA")
group<-read.csv('../08_GSEA/risk.csv')
colnames(group) <- c('sample', 'group')
group<-group[order(group$group),]

ssgsea_score<-read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1)
ssgsea_score<-ssgsea_score[,group$sample]
annotation_col<-as.data.frame(group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)
ssgsea_score <- as.matrix(ssgsea_score)

library(pheatmap)
color.key<-c("#69bcce","#fed53b")
ann_colors<-list(Group = c('High Risk'="#69bcce",'Low Risk'="#fed53b"))
p <- pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
p
png(filename = "01.heatmap.png", height = 5, width = 8, units = 'in',res = 600,family='Times')
p
dev.off()
pdf(file = "01.heatmap.pdf", height = 5,width = 8,family='Times')
p
dev.off()



rm(list = ls())

setwd("./09_ssGSEA")
group<-read.csv('../08_GSEA/risk.csv')
colnames(group) <- c('sample', 'group')
group_case<-group[group$group=='High Risk',]
group_case<-as.character(group_case$sample)
group_control<-group[group$group=='Low Risk',]
group_control<-as.character(group_control$sample)
library(dplyr)
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) %>% lc.tableToNum
tiics_result <- tiics_result[, group$sample] %>% as.matrix()

pvalue = padj <- matrix(0, nrow(tiics_result), 1)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value


}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(
  pvalue, 
  padj,
  row.names = rownames(tiics_result))









rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$pvalue<0.05),]
dim(diff_Table)

write.csv(rTable,
          file = "02.cell_all_wilcox_test.csv",
          quote = F,
          row.names = F)
write.csv(diff_Table,
          file = "02.cell_diff_wilcox_test.csv",
          quote = F,
          row.names = F)
diff_Table$immune_cell






library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

cell.data <- data.frame(Immune_Cell=rownames(tiics_result), 
                        tiics_result, 
                        pvalue=rTable$pvalue)
plot.cell <- cell.data[which(cell.data$pvalue<0.05),]

diff_tiics <- rownames(plot.cell)
violin_dat <- gather(plot.cell, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_control,
                           "Low Risk", "High Risk") 
violin_dat$Group <- factor(violin_dat$Group, levels = c("Low Risk", "High Risk"))
violin_dat <- violin_dat[,-2]
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  geom_violin(trim=T,color='black') + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,
               color='black')+ 



  scale_fill_manual(values= c("#69bcce90","#fed53b90"))+ #设置填充的颜色
  labs(title="Cell Infiltration", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     method = "wilcox.test", 
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), 
        axis.text.x=element_text(angle=50,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=14,face="bold"), 
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), 
        legend.title=element_text(face="bold", colour="black", size=11),
        text=element_text(family = 'Times'),


        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
boxplot_diff_TIICs
ggsave('02.ssgsea_Box.pdf',boxplot_diff_TIICs,w=8,h=5)
ggsave('02.ssgsea_Box.png',boxplot_diff_TIICs,w=8,h=5)







rm(list = ls())

setwd("./09_ssGSEA")












res3 <- read.csv('01.ssgsea_result_cell.csv',row.names = 1)
dif_cell <- read.csv('02.cell_diff_wilcox_test.csv')
res3 <- res3[rownames(res3) %in% dif_cell$immune_cell,]

res3 <- t(res3) %>% as.data.frame()






library(ggcorrplot)
cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)

dat.r <- data.frame(cor_data)
dat.r$cell <- rownames(dat.r)
dat.r <- gather(dat.r,'Gene','r',-c(cell))

dat.p <- data.frame(corp)
dat.p$cell <- rownames(dat.p)
dat.p <- gather(dat.p,'Gene','p',-c(cell))
corrp <- (cbind(dat.r, dat.p))[,c("cell","Gene",'r',"p")]
colnames(corrp) <- c("cell1","cell2",'Correlation',"Pvalue")

corrp$cell2 <- gsub("\\.", " ",corrp$cell2)
write.csv(corrp,file="03.cell_cor.csv",quote=F)









library(ggcorrplot)
library(corrplot)
pdf("03.cell_corHeatmap.pdf",height=14,width=14,family='Times')
col1 <- colorRampPalette(c("#69bcce", "white", "#fed53b"))
corr_plot<-corrplot(cor_data,
                    tl.cex=1.5,
                    method = "square", 
                    is.corr = T,
                    type = "lower", 
                    p.mat = corp,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 0.8,
                    number.cex = 1.2
)
dev.off()

png("03.cell_corHeatmap.png",height=14,width=14,family='Times',units='in',res=600)
col1 <- colorRampPalette(c("#69bcce", "white", "#fed53b"))
corr_plot<-corrplot(cor_data,
                    tl.cex=1.5,
                    method = "square", 
                    is.corr = T,
                    type = "lower", 
                    p.mat = corp,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 0.8,
                    number.cex = 1.2
)
dev.off()

























































































































































































































































































































