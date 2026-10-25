
rm(list = ls())

setwd("./10_TIDE")
tcga_mrna<-read.csv("../00_rawdata/output_mRNA_lncRNA_expr/TCGA-HNSC_mrna_expr_tpm.csv",check.names = F,row.names = 1)
new_colnames <- substr(colnames(tcga_mrna), 1, 15) 
colnames(tcga_mrna) <- new_colnames  

dat_score<- read.csv('../05_model/01.OSCC_train_risk.csv') %>% lc.tableToNum()
fpkm <- tcga_mrna[,dat_score$X]






rownmean <- apply(fpkm,1,mean)
data<- sweep(fpkm, 1, rownmean)
write.table(data,'TIDE_rawdata.xls',sep='\t',quote=F)  


tide<-read.csv('TIDE.csv')
tide<-data.frame(tide[,c(1,4)])
colnames(tide)[1] <- 'sample'

group<-read.csv('../08_GSEA/risk.csv')
colnames(group) <- c('sample', 'risk')







dat <- merge(tide, group, by= 'sample')
write.csv(dat, file = '01.TIDEscore.csv', row.names = F, quote = F)

data <- na.omit(dat)
label.y = c(2.2,0.7,2.2)
plot.ips <- c('TIDE', 'Dysfunction', 'Exclusion')
for (i in c(1)) {

  my_comparisons = list( c("High Risk", "Low Risk"))
  p <- ggplot(data, aes(risk, get(plot.ips[i]))) + 
    ggdist::stat_halfeye(aes(color=risk,fill=risk),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(aes(color=risk),width = .1, outlier.shape = NA) +



    geom_jitter(aes(color=risk),width = .05, alpha = .3) +
    scale_fill_manual(values=c("#f4c466","#7974a5")) +
    scale_color_manual(values=c("#f4c466","#7974a5")) + 
    coord_flip()+
    theme_bw()+
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",label.y = label.y[i],
                       color='black',family="Times",face = "bold")+
    labs(x = "risk", y = plot.ips[i], title = "") + 
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
          axis.text.y = element_text(size = 0, face = "bold", family = "Times",color='black'),
          legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
          legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
          legend.position = "bottom"
    )+
    theme(panel.grid =element_blank()) +   
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_blank())
  
  png(paste0('02.',plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0('02.',plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}

rm(list = ls())

setwd("./10_TIDE")
tide<-read.csv('TIDE.csv')
colnames(tide)
tide<-data.frame(tide[,c("Patient","Dysfunction")])

colnames(tide)[1] <- 'sample'

group<-read.csv('../08_GSEA/risk.csv')
colnames(group) <- c('sample', 'risk')

dat <- merge(tide, group, by= 'sample')
write.csv(dat, file = '02.Dysfunction.csv', row.names = F, quote = F)

data <- na.omit(dat)
label.y = c(2.2,0.7,2.2)
plot.ips <- c('Dysfunction')
for (i in c(1)) {

  my_comparisons = list( c("High Risk", "Low Risk"))
  p <- ggplot(data, aes(risk, get(plot.ips[i]))) + 
    ggdist::stat_halfeye(aes(color=risk,fill=risk),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(aes(color=risk),width = .1, outlier.shape = NA) +



    geom_jitter(aes(color=risk),width = .05, alpha = .3) +
    scale_fill_manual(values=c("#f4c466","#7974a5")) +
    scale_color_manual(values=c("#f4c466","#7974a5")) + 
    coord_flip()+
    theme_bw()+
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",label.y = label.y[i],
                       color='black',family="Times",face = "bold")+
    labs(x = "risk", y = plot.ips[i], title = "") + 
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
          axis.text.y = element_text(size = 0, face = "bold", family = "Times",color='black'),
          legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
          legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
          legend.position = "bottom"
    )+
    theme(panel.grid =element_blank()) +   
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_blank())
  
  png(paste0('03.',plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0('03.',plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}


rm(list = ls())

setwd("./10_TIDE")
tide<-read.csv('TIDE.csv')
colnames(tide)
tide<-data.frame(tide[,c("Patient","Exclusion")])

colnames(tide)[1] <- 'sample'

group<-read.csv('../08_GSEA/risk.csv')
colnames(group) <- c('sample', 'risk')

dat <- merge(tide, group, by= 'sample')
write.csv(dat, file = '02.Exclusion.csv', row.names = F, quote = F)

data <- na.omit(dat)
label.y = c(2.2,0.7,2.2)
plot.ips <- c('Exclusion')
for (i in c(1)) {

  my_comparisons = list( c("High Risk", "Low Risk"))
  p <- ggplot(data, aes(risk, get(plot.ips[i]))) + 
    ggdist::stat_halfeye(aes(color=risk,fill=risk),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(aes(color=risk),width = .1, outlier.shape = NA) +



    geom_jitter(aes(color=risk),width = .05, alpha = .3) +
    scale_fill_manual(values=c("#f4c466","#7974a5")) +
    scale_color_manual(values=c("#f4c466","#7974a5")) + 
    coord_flip()+
    theme_bw()+
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",label.y = label.y[i],
                       color='black',family="Times",face = "bold")+
    labs(x = "risk", y = plot.ips[i], title = "") + 
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
          axis.text.y = element_text(size = 0, face = "bold", family = "Times",color='black'),
          legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
          legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
          legend.position = "bottom"
    )+
    theme(panel.grid =element_blank()) +   
    theme(axis.text = element_blank()) +  
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_blank())
  
  png(paste0('04.',plot.ips[i],'.png'),width=7,height=6,family='Times',units='in',res=600)
  print(p)
  dev.off()
  pdf(paste0('04.',plot.ips[i],'.pdf'),width=7,height=6,family='Times')
  print(p)
  dev.off()
}
