rm(list = ls())

setwd('./06_progmodel')



survival <- read.csv('../02_uncox_pan/OSCC_OS.csv',check.names = F)
survival <- survival %>% na.omit()
colnames(survival)<-c("sample","OS","OS.time")
phenotype <- read_tsv(file = '../00_rawdata/TCGA-HNSC.GDC_phenotype.tsv.gz')


train_phenotype <- data.frame(sample = phenotype$submitter_id.samples,
                              age = phenotype$age_at_initial_pathologic_diagnosis,
                              gender = phenotype$gender.demographic,
                              race = phenotype$race.demographic,
                              TNM.stage = phenotype$tumor_stage.diagnoses,



                              grade = phenotype$neoplasm_histologic_grade)


train_phenotype$sample <- substr(train_phenotype$sample, 1, nchar(train_phenotype$sample) - 1)

train_phenotype <- merge(train_phenotype, survival, by='sample')

write.csv(train_phenotype, file = '01.phenotype.csv', row.names = F, quote = F)



train_phenotype2 <- train_phenotype
train_phenotype2$OS <- as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time <- as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$age) 
train_phenotype2$age <- as.numeric(train_phenotype2$age)


table(train_phenotype2$gender)
train_phenotype2$gender <- ifelse(train_phenotype2$gender == 'female',1,0)

table(train_phenotype2$race)
train_phenotype2$race <- gsub('american indian or alaska native',NA,train_phenotype2$race)
train_phenotype2$race <- gsub('asian','1',train_phenotype2$race)
train_phenotype2$race <- gsub('black or african american','2',train_phenotype2$race)
train_phenotype2$race <- gsub('white','3',train_phenotype2$race)
train_phenotype2$race <- gsub('not reported',NA,train_phenotype2$race)

table(train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('a','',train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('b','',train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('c','',train_phenotype2$TNM.stage)

train_phenotype2$TNM.stage <- gsub('not reported',NA,train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('stge iii','3',train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('stge iv','4',train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('stge ii','2',train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage <- gsub('stge i','1',train_phenotype2$TNM.stage)
































table(train_phenotype2$grade)
train_phenotype2$grade <- gsub('G3','3',train_phenotype2$grade)
train_phenotype2$grade <- gsub('G4',NA,train_phenotype2$grade)
train_phenotype2$grade <- gsub('G2','2',train_phenotype2$grade)
train_phenotype2$grade <- gsub('G1','1',train_phenotype2$grade)
train_phenotype2$grade <- gsub('GX',NA,train_phenotype2$grade)


riskScore <- read.csv('../05_model/01.OSCC_train_risk.csv') %>% lc.tableToNum()
colnames(riskScore)
riskScore <- subset(riskScore, select = c(X, riskScore))
colnames(riskScore)<-c("id","riskScore")

colnames(train_phenotype2)
colnames(train_phenotype2)[1] <- 'id'

os_risk_clinical <- merge(train_phenotype2, riskScore, by = "id") %>% 
  column_to_rownames(var = 'id') %>% lc.tableToNum()
dim(os_risk_clinical)


colnames_os_risk_clinical <- colnames(os_risk_clinical)
covariates <- colnames_os_risk_clinical[-which(colnames_os_risk_clinical %in% c("OS", "OS.time"))]


os_risk_clinical$gender <- factor(os_risk_clinical$gender)
os_risk_clinical$race <- factor(os_risk_clinical$race)
os_risk_clinical$TNM.stage <- factor(os_risk_clinical$TNM.stage)



os_risk_clinical$grade <- factor(os_risk_clinical$grade)




library(survival)
library(survminer)
res.all <- data.frame()
for (i in 1:length(covariates)) {

  gene_list <- covariates[i]
  cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                                paste(gene_list,
                                      sep = '',
                                      collapse = '+')))
  
  cox_more <- coxph(cox_data,
                    data = os_risk_clinical)
  cox_zph <- cox.zph(cox_more)
  cox_table <- cox_zph$table[-nrow(cox_zph$table),]
  cox_table
  res <- data.frame(gene=covariates[i],PH.p=cox_table[3])
  res.all <- rbind(res.all, res)
  p <- ggcoxzph(cox_zph)







}
colnames(res.all)[1] <- 'phenotype'
write.csv(res.all,
          file = "01.uniCox_PH.csv",
          quote = F,
          row.names = F)

ph.gene <- res.all[as.numeric(res.all$PH.p)>0.05,]
dim(ph.gene)
ph.gene$phenotype





library(survival)

res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = os_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])

res.age = coxph(Surv(time = OS.time, event = OS) ~ age, data = os_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ gender, data = os_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

res.race = coxph(Surv(time = OS.time, event = OS) ~ race, data = os_risk_clinical) %>% summary
res.race = cbind(res.race$conf.int[,-2], res.race$coefficients[,5])

res.TNM.stage = coxph(Surv(time = OS.time, event = OS) ~ TNM.stage, data = os_risk_clinical) %>% summary
res.TNM.stage = cbind(res.TNM.stage$conf.int[,-2], res.TNM.stage$coefficients[,5])










res.grade = coxph(Surv(time = OS.time, event = OS) ~ grade, data = os_risk_clinical) %>% summary
res.grade = cbind(res.grade$conf.int[,-2], res.grade$coefficients[,5])


res.ref = c(1,1,1,NA)

res = rbind(res.risk, res.age,res.gender,  res.ref, res.TNM.stage, res.ref, res.grade) %>% as.data.frame()
rownames(res)
res$Indicators = c('riskScore',
                   'age',
                   'gender', 

                   'Stage i(Reference)','Stage ii','Stage iii','Stage iv',



                   'G1(Reference)', 'G2', 'G3')

colnames(res) = c("hr","low","up","pv","Indicator")




res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
res2 <- subset(res2, select = -c(Indicator))
write.csv(res2, file = "02.univariate_cox.csv", quote = F, row.names = T)



library(survival)
library(survminer)
colnames(os_risk_clinical)
cox_gene <- c("riskScore","age","TNM.stage","gender","grade")
a<-c("OS","OS.time",cox_gene)
train_dat_cox_gene_sv<-os_risk_clinical[,colnames(os_risk_clinical)%in%a]


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

  p <- ggcoxzph(cox_zph)
  pdf(paste0('01.ph.test_',gene_list,'.pdf'), width = 6,height = 4,onefile = F,family='Times')
  print(p)
  dev.off()
  
  png(paste0('01.ph.test_',gene_list,'.png'), width = 6,height = 4,family='Times',units='in',res=600)
  print(p)
  dev.off()
}


library(tidyr)

hz <- paste(round(res2$HR,3), 
            "(",
            round(res2$HR.95L,3),
            "-",
            round(res2$HR.95H,3),
            ")",
            sep = "") 
hz

hz[c(4,8)] <- "" 

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.05,"< 0.05",round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1


library(forestplot)
pdf(file = '02.uncoxforest.pdf',height = 6, width = 12, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,rep(FALSE, 70)),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), 
           upper=c(NA,res2$HR.95H), 
           boxsize=0.15,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=2,    
           colgap=unit(5,"mm"),    
           xticks = c(0.5,1,15,18), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate Cox") 
dev.off()

png(file = '02.uncoxforest.png',height = 6, width = 12, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,rep(FALSE, 70)),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), 
           upper=c(NA,res2$HR.95H), 
           boxsize=0.15,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=2,    
           colgap=unit(5,"mm"),    
           xticks = c(0.5,1,15,18), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate Cox") 
dev.off()


cox_more <- coxph(Surv(time = OS.time, event = OS) ~ riskScore+TNM.stage+age, data = os_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]      
cox_table





res.mul <- coxph(Surv(time = OS.time, event = OS) ~ riskScore+TNM.stage+age, data = os_risk_clinical) %>% summary
res.mul  <-  cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)


res.mul$Indicators = c("riskScore",'Stage ii','Stage iii','Stage iv','age')
colnames(res.mul) = c("hr","low","up","pv","Indicator")





multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
multi_res <- subset(multi_res, select = -c(Indicator))
write.csv(multi_res,
          file = "03.multivariate_cox.csv",
          quote = F,
          row.names = T)

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[2] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.05,
                                      "< 0.05",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '03.mulcoxforest.pdf',height = 6, width = 12, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), 
           upper=c(NA,multi_res$HR.95H), 
           boxsize=0.1,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),    
           xticks = c(0.5,1,5,10,15,19), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate Cox") 
dev.off()

png(file = '03.mulcoxforest.png',height = 6, width = 12, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=2,  
           is.summary = c(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="#d9b02c", lines="#4981bb", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), 
           upper=c(NA,multi_res$HR.95H), 
           boxsize=0.1,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),    
           xticks = c(0.5,1,5,10,15,19), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawNormalCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),





           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate Cox") 
dev.off()


































































































































































































































































































































































