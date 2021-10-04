#!/usr/bin/env Rscript

#Trait transformation and stratify data by sex
library(data.table)

## INPUT
args = commandArgs(trailingOnly=TRUE)
cat(args, sep = "\n")

Combined <- fread(args[1])
MALES <- subset(Combined, Sex=="1")
FEMALES <- subset(Combined, Sex=="0")

#linear model
MALESlm <-lm(BMI ~ Age + Disease_status, data=MALES)
FEMALESlm <-lm(BMI ~ Age + Disease_status, data=FEMALES)
Combinedlm <-  lm(BMI ~ Age + Sex + Disease_status, data=Combined)


#adding residuals to the cov file
MALES$residuals=MALESlm$residuals
FEMALES$residuals=FEMALESlm$residuals
Combined$residuals=Combinedlm$residuals

#inverse rank normalisation of traits
MALES$pheno <- qnorm((rank(MALES$residuals,na.last="keep") - 0.5)/sum(!is.na(MALES$residuals)))
FEMALES$pheno <- qnorm((rank(FEMALES$residuals,na.last="keep") - 0.5)/sum(!is.na(FEMALES$residuals)))
Combined$pheno <- qnorm((rank(Combined$residuals,na.last="keep") - 0.5)/sum(!is.na(Combined$residuals)))

#Compute variance of BMI
MALES$varpheno <- (MALES$pheno)^2
FEMALES$varpheno <- (FEMALES$pheno)^2
Combined$varpheno <- (Combined$pheno)^2

#create working pheno file
FEMALES.pheno <- subset(FEMALES, select=c("FID","IID","pheno","varpheno"))
MALES.pheno <- subset(MALES,select=c("FID","IID","pheno","varpheno"))
Combined.pheno <- subset(Combined,select=c("FID","IID","pheno","varpheno"))

write.table(MALES.pheno, file = "MALES.pheno", quote =F,col.names =T,row.names=F,sep = " ")
write.table(FEMALES.pheno, file = "FEMALES.pheno", quote =F,col.names =T,row.names=F,sep = " ")
write.table(Combined.pheno, file = "Combined.pheno", quote =F,col.names =T,row.names=F,sep = " ")
