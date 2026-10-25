rm(list = ls())

setwd("./00_rawdata")
library(easyTCGA)

tcga_expr_file <- "tcga_RSEM_gene_tpm.gz"
tcga_clin_file <- "Survival_SupplementalTable_S1_20171025_xena_sp"



getpancancer_xena(tcga_expr_file = tcga_expr_file,
                  tcga_clin_file = tcga_clin_file,
                  type = "tcga")


load(file="output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")



tcga_mrna <- tcga_mrna_clin[, c(1, 3, 35:ncol(tcga_mrna_clin))]

save(tcga_mrna, file = "tcga_mrna.rdata")


print("tcga_mrna has been saved to tcga_mrna.rdata")



head(tcga_mrna[, 1:4], 4)


CMGRs <- read.csv('CMGRs.csv',stringsAsFactors = FALSE)
genes_of_interest<-CMGRs$SYMBOL

selected_columns <- c("sample_id", "project", genes_of_interest)

tcga_mrna_clin_subset <- tcga_mrna_clin[, colnames(tcga_mrna_clin) %in% selected_columns]
























missing_columns <- setdiff(selected_columns, colnames(tcga_mrna_clin_subset))

print(missing_columns)




head(tcga_mrna_clin[, 1:4], 4)
write.csv(tcga_mrna_clin_subset,file = 'tcga.tpm.csv')
colnames(tcga_mrna_clin)
phenotype <- cbind(tcga_mrna_clin[,c(1:34)])
phenotype <- phenotype[,c(1:3,26:33)]
write.csv(phenotype,file = 'survival.tcga.csv')


phe1 <- read.delim2("../00_rawdata/TCGA-HNSC.GDC_phenotype.tsv.gz", sep = "\t")
phe2 <- subset(phe1,site_of_resection_or_biopsy.diagnoses=='Anterior floor of mouth'|
                 site_of_resection_or_biopsy.diagnoses=='Base of tongue, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Border of tongue'|
                 site_of_resection_or_biopsy.diagnoses=='Cheek mucosa'|
                 site_of_resection_or_biopsy.diagnoses=='Floor of mouth, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Gum, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Hard palate'|
                 site_of_resection_or_biopsy.diagnoses=='Lip, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Lower gum'|
                 site_of_resection_or_biopsy.diagnoses=='Mandible'|
                 site_of_resection_or_biopsy.diagnoses=='Mouth, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Palate, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Tongue, NOS'|
                 site_of_resection_or_biopsy.diagnoses=='Upper gum'|
                 site_of_resection_or_biopsy.diagnoses=='Ventral surface of tongue, NOS')
id<-phe2$submitter_id.samples
phenotype <- cbind(tcga_mrna_clin[,c(1:34)])
phenotype <- phenotype[,c(1:3,26:33)]




id <- substr(phe2$submitter_id.samples, 1, nchar(phe2$submitter_id.samples) - 1)

OSCC_survival <- phenotype[phenotype$sample_id%in% id, ]
OSCC_tpm<- tcga_mrna_clin_subset[tcga_mrna_clin_subset$sample_id%in% id, ]
tumor_type <- ifelse(substr(OSCC_tpm$sample_id, 14, 14) == "0", "肿瘤样本", "正常样本")
table(tumor_type)
OSCC_survival <- phenotype[phenotype$sample_id%in% id, ]


write.csv(id,file = 'OSCC.id.csv', row.names = FALSE)

write.csv(OSCC_survival,file = 'survival.OSCC.csv', row.names = FALSE)
OSCC_tpm<- tcga_mrna_clin_subset[tcga_mrna_clin_subset$sample_id%in% id, ]
write.csv(OSCC_tpm,file = 'OSCC.tpm.csv', row.names = FALSE)





getmrnaexpr("TCGA-HNSC")
dat_OSCC1<-read.csv("output_mRNA_lncRNA_expr/TCGA-HNSC_mrna_expr_count.csv",check.names = F)

OSCC1 <- dat_OSCC1[colnames(dat_OSCC1)%in% id, ]














rm(list = ls())

setwd('./00_rawdata')

library(GEOquery)
library(Biobase)
gset<-getGEO('GSE41613',
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gset

gpl<-getGEO("GPL570",destdir = '.')

gpl<-Table(gpl)
colnames(gpl)



gpl$`Gene Symbol`
library(tidyverse)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter(`Gene Symbol`!='')%>%
  separate('Gene Symbol',c('Symbol','drop'),sep = ' /// ')%>%
  dplyr::select(-drop)


probe2symobl=probe2symobl[probe2symobl$Symbol!='',] %>% na.omit()
probe2symobl <- probe2symobl[!probe2symobl$Symbol == "previous version conserved probe", ]


dat_expr<-expr
dat_expr$ID<-rownames(dat_expr)
dat_expr$ID<-as.character(dat_expr$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat_expr<-dat_expr %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(Symbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(Symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   

pcg <- read.delim2('/data/nas3/mayanyu/project/01_YQZZ-10417-4/PCG(1).xls(v22)')

dat_expr <- dat_expr[rownames(dat_expr)%in%pcg$gene_name,] %>% na.omit()


dim(dat_expr)

pd<-pData(gset[[1]])
colnames(pd)

group <- data_frame(sample = rownames(pd), OS= pd$characteristics_ch1.4, OS.time= pd$characteristics_ch1.5)







group1 <- group %>%
  separate(OS, into = c("type1", "patient_group"), sep = "-", extra = "merge") %>%
  dplyr::select(sample, type1, patient_group) %>%
  rename(OS = type1)
group1$OS <- ifelse(group1$OS=='vital: Alive', '0', '1')
group$OS.time <- str_extract(group$OS.time, "\\d+\\.?\\d*")
group$OS.time <- as.numeric(group$OS.time)
group$OS.time <- group$OS.time * 365
group$OS<-group1$OS 


























identical(group$sample, colnames(dat_expr))


write.csv(dat_expr, file = 'dat(GSE41613).csv',row.names = T)
write.csv(group, file = 'OS(GSE41613).csv',row.names = F)













