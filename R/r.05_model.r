rm(list = ls())

setwd("./05_model")


library(lance)
library(readxl)
library(readr)
library(tidyverse)
rm(list = ls())
data0 <- read.csv('../02_uncox_pan/OSCC_TPM.csv')
survival <- read.csv('../02_uncox_pan/OSCC_OS.csv')
gene<- read.csv("../04_OSCC_uniCox/01.OSCC_gene.csv")
matched_indices <- c("sample_id",gene$symbol)  
data1<-data0[,colnames(data0)%in%matched_indices]
train_dat<-merge(survival,data1,by='sample_id')
rownames(train_dat)<-train_dat$sample_id
train_dat<-train_dat[,-1]
colnames(train_dat)
train_dat<-na.omit(train_dat)
train_dat <- train_dat[train_dat$OS.time > 0, ]


survival_test<-read.csv("../00_rawdata/OS(GSE41613).csv",header=TRUE)
survival_test <- survival_test[survival_test$OS.time>0,]
data_test <- read.csv('../00_rawdata/dat(GSE41613).csv',header=TRUE, row.names=1)
survival_dat<-t(data_test)
survival_dat<-survival_dat[rownames(survival_dat)%in%survival_test$sample,] 
test_dat<-survival_dat[,colnames(survival_dat)%in%matched_indices]
test_dat<-as.data.frame(test_dat)
test_dat$sample<-rownames(test_dat)
test_dat<-merge(survival_test,test_dat,by='sample')
rownames(test_dat)<-test_dat$sample
test_dat<-test_dat[,-1]
test_dat <- na.omit(test_dat)
colnames(test_dat)



x_train <- subset(train_dat, select = -c(OS, OS.time))
y_train <- subset(train_dat, select = c(OS, OS.time))
x_test <- subset(test_dat, select = -c(OS, OS.time))
y_test <- subset(test_dat, select = c(OS, OS.time))

library(randomForestSRC)
library(e1071)
library(survival)
library(caret)
library(survivalROC)
library(ggplot2)


set.seed(1)


x_train <- scale(x_train)
x_test <- scale(x_test)


x_test <- x_test[, colnames(x_train)]


train_data <- data.frame(y_train, x_train)
test_data <- data.frame(y_test, x_test)











































































set.seed(1)  
final_model <- rfsrc(Surv(OS.time, OS) ~ ., data = train_data,
                     ntree = 100, mtry = 4)


train_predictions <- predict(final_model, newdata = train_data)$predicted
test_predictions <- predict(final_model, newdata = test_data)$predicted


train_data$riskScore <- train_predictions
test_data$riskScore <- test_predictions
write.csv(train_data,file = '01.OSCC_train_risk.csv',row.names = T)
write.csv(test_data,file = '01.OSCC_test_risk.csv',row.names = T)



rm(list = ls())


setwd("./05_model")

train_dat<-read.csv("01.OSCC_train_risk.csv",row.names = 1,check.names = F)

dat_score <- train_dat[, c('OS', 'OS.time', 'riskScore')]
set.seed(1) 
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS',variables = 'riskScore')
cutpoint <- res.cut$cutpoint[1]
cutpoint
dat_score$risk <- ifelse(dat_score$riskScore > 39.69651, 1, 0)
dat_score$risk <- ifelse(dat_score$risk == 1, 'high risk', 'low risk')
table(dat_score$risk)

kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = dat_score)
library(survival)
library(survminer)


km_plot <- ggsurvplot(kmfit,
                      pval = TRUE,
                      conf.int = FALSE,
                      legend.labs = c("High risk", "Low risk"),
                      legend.title = "Risk score",
                      title = "Train Data KM Curve",
                      font.main = c(15, "bold"),
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      linetype = "strata",
                      ggtheme = theme_bw(),
                      palette = c("#f47e84","#6694e9"))


combined_plot <- ggarrange(km_plot$plot, km_plot$table, ncol = 1, heights = c(3, 1))


ggsave("01.Train_KM.pdf", plot = combined_plot, width = 8, height = 8)


ggsave("01.Train_KM.png", plot = combined_plot, width = 8, height = 8, dpi = 600)
library(timeROC)
library(survival)


df<-dat_score
ROC <- timeROC(T=df$OS.time,   
               delta=df$OS,   
               marker=df$riskScore,   
               cause=1,                
               weighting="marginal",   
               times=c(1*365, 3*365, 5*365),       
               iid=TRUE)
ROC
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)


legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],2)), 
         paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n") 
pdf(file = "02.Train_ROC_Curve.pdf", width = 8, height = 8)
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 1 year: ", round(ROC[["AUC"]][1], 2)), 
         paste0("AUC at 3 year: ", round(ROC[["AUC"]][2], 2)), 
         paste0("AUC at 5 year: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n") 
dev.off()


png(file = "02.Train_ROC_Curve.png", width = 600, height = 600)
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 1 year: ", round(ROC[["AUC"]][1], 2)), 
         paste0("AUC at 3 year: ", round(ROC[["AUC"]][2], 2)), 
         paste0("AUC at 5 year: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n") 
dev.off()





risk<-dat_score
id<-rownames(dat_score)

risk$risk <- factor(risk$risk, levels = c("high risk", "low risk"), labels = c("High Risk", "Low Risk"))
library(ggplot2)

risk_dis <- ggplot(risk, aes(x = reorder(id, riskScore),
                             y = riskScore,
                             color = factor(risk,
                                            levels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("High Risk" = "#f47e84", "Low Risk" = "#6694e9")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)],
                   labels = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  geom_hline(yintercept = cutpoint$cutpoint,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Risk Score",
       title = "Train") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  
        legend.background = element_rect(color = "black", linewidth = 0.3),  
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(risk_dis)


surv_stat <- ggplot(risk, aes(x = reorder(id, riskScore),
                              y = OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point(na.rm = F) +  
  scale_color_manual(values = c("#f47e84","#6694e9")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)],
                   labels = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100),
                   expand = c(0.02, 0)) +
  ylim(c(0, max(risk$OS.time / 365, na.rm = TRUE))) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Survival time (Years)",
       title = "") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  
        legend.background = element_rect(color = "black", linewidth = 0.3),  
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(surv_stat)


library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)

ggsave(filename = '03.Train_risk_and_survival_Distribution.pdf',p,w=6,h=7)
ggsave(filename = '03.Train_risk_and_survival_Distribution.png',p,w=6,h=7, units='in', dpi = 600)

library(ggplot2)


cutoff_value <- 39.69651


density_with_histogram <- ggplot(risk, aes(x = riskScore)) +

  geom_density(fill = "#69b3a2", color = "#1f77b4", alpha = 0.6) +

  geom_histogram(aes(y = ..density..), 
                 binwidth = 2, 
                 fill = "
                 color = "
                 alpha = 0.4) +

  geom_vline(xintercept = cutoff_value, 
             linetype = "dashed", 
             color = "black", 
             size = 1) +

  labs(x = "Risk Score",
       y = "Density"

  )
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(size = 12),
        text = element_text(size = 10))


print(density_with_histogram)



ggsave(filename = "04.Train_density_with_histogram.pdf", 
       plot = density_with_histogram, 
       width = 8, 
       height = 6, 
       units = "in")


ggsave(filename = "04.Train_density_with_histogram.png", 
       plot = density_with_histogram, 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)

















rm(list = ls())


setwd("./05_model")
train_dat<-read.csv("01.OSCC_test_risk.csv",row.names = 1,check.names = F)

dat_score <- train_dat[, c('OS', 'OS.time', 'riskScore')]
set.seed(1) 
res.cut <- surv_cutpoint(dat_score, time = 'OS.time', event = 'OS',variables = 'riskScore')
cutpoint <- res.cut$cutpoint[1]
cutpoint
dat_score$risk <- ifelse(dat_score$riskScore > 44.80562, 1, 0)
dat_score$risk <- ifelse(dat_score$risk == 1, 'high risk', 'low risk')
table(dat_score$risk)

kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = dat_score)
library(survival)
library(survminer)


km_plot <- ggsurvplot(kmfit,
                      pval = TRUE,
                      conf.int = FALSE,
                      legend.labs = c("High risk","Low risk"),
                      legend.title = "Risk score",
                      title = " ",
                      font.main = c(15, "bold"),
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      linetype = "strata",
                      ggtheme = theme_bw(),
                      palette = c("#f47e84","#6694e9"))


combined_plot <- ggarrange(km_plot$plot, km_plot$table, ncol = 1, heights = c(3, 1))


ggsave("05.Test_KM.pdf", plot = combined_plot, width = 8, height = 8)


ggsave("05.Test_KM.png", plot = combined_plot, width = 8, height = 8, dpi = 600)
library(timeROC)
library(survival)


df<-dat_score
ROC <- timeROC(T=df$OS.time,   
               delta=df$OS,   
               marker=df$riskScore,   
               cause=1,                
               weighting="marginal",   
               times=c(1*365, 3*365, 5*365),       
               iid=TRUE)
ROC
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)


legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],2)), 
         paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n") 
pdf(file = "06.Test_ROC_Curve.pdf", width = 8, height = 8)
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 1 year: ", round(ROC[["AUC"]][1], 2)), 
         paste0("AUC at 3 year: ", round(ROC[["AUC"]][2], 2)), 
         paste0("AUC at 5 year: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n") 
dev.off()


png(file = "06.Test_ROC_Curve.png", width = 600, height = 600)
plot(ROC, 
     time=1*365, col="red", lwd=2, title = "")   
plot(ROC,
     time=3*365, col="blue", add=TRUE, lwd=2)    
plot(ROC,
     time=5*365, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 1 year: ", round(ROC[["AUC"]][1], 2)), 
         paste0("AUC at 3 year: ", round(ROC[["AUC"]][2], 2)), 
         paste0("AUC at 5 year: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2, bty="n") 
dev.off()





risk<-dat_score
id<-rownames(dat_score)

risk$risk <- factor(risk$risk, levels = c("high risk", "low risk"), labels = c("High Risk", "Low Risk"))
library(ggplot2)

risk_dis <- ggplot(risk, aes(x = reorder(id, riskScore),
                             y = riskScore,
                             color = factor(risk,
                                            levels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("High Risk" = "#f47e84", "Low Risk" = "#6694e9")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)],
                   labels = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100),
                   expand = c(0.02, 0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  geom_hline(yintercept = cutpoint$cutpoint,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Risk Score",
       title = " ") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  
        legend.background = element_rect(color = "black", linewidth = 0.3),  
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(risk_dis)


surv_stat <- ggplot(risk, aes(x = reorder(id, riskScore),
                              y = OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point(na.rm = F) +  
  scale_color_manual(values = c("#6694e9", "#f47e84")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100)],
                   labels = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100),
                   expand = c(0.02, 0)) +
  ylim(c(0, max(risk$OS.time / 365, na.rm = TRUE))) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == "Low Risk"),]) + 0.5,
             linetype = "dashed") +
  labs(x = "Patients (increasing risk score)",
       y = "Survival time (Years)",
       title = "") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",  
        legend.background = element_rect(color = "black", linewidth = 0.3),  
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(family = 'Times'),
        text = element_text(family = 'Times'))

print(surv_stat)


library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)
ggsave(filename = '07.Test_risk_and_survival_Distribution.pdf',p,w=6,h=7)
ggsave(filename = '07.Test_risk_and_survival_Distribution.png',p,w=6,h=7, units='in', dpi = 600)

library(ggplot2)


cutoff_value <- 44.80562


density_with_histogram <- ggplot(risk, aes(x = riskScore)) +

  geom_density(fill = "#69b3a2", color = "#1f77b4", alpha = 0.6) +

  geom_histogram(aes(y = ..density..), 
                 binwidth = 2, 
                 fill = "
                 color = "
                 alpha = 0.4) +

  geom_vline(xintercept = cutoff_value, 
             linetype = "dashed", 
             color = "black", 
             size = 1) +

  labs(x = "Risk Score",
       y = "Density"

  )
theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(size = 12),
        text = element_text(size = 10))


print(density_with_histogram)



ggsave(filename = "08.Test_density_with_histogram.pdf", 
       plot = density_with_histogram, 
       width = 8, 
       height = 6, 
       units = "in")


ggsave(filename = "08.Test_density_with_histogram.png", 
       plot = density_with_histogram, 
       width = 8, 
       height = 6, 
       units = "in", 
       dpi = 300)




