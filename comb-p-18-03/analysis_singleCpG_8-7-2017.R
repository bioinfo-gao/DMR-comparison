setwd("C:/Users/lxw391/Dropbox (BBSR)/Zhen_Gao/LW-test/AD-combp")

# AD analysis - get single CpG p-values for comb-p 

# beta values
library(DMRcate)
betas <- readRDS("C:/Users/lxw391/Dropbox (BBSR)/Zhen_Gao/LW_Projects/For-Lily/2017-07-correlation-and-map/Annot/betaValue.rds")   
  
mvalues <- log2(betas/(1-betas))
mvalues.filtered <- rmSNPandCH(as.matrix(mvalues), dist=2, mafcut=0.05, rmXY=TRUE, rmcrosshyb = TRUE, and=FALSE)

# pheno info
library(openxlsx)
phenofile <- "C:/Users/lxw391/Box Sync/METHODS-SHARED/METHOD-METHYL-GENE-TEST/AD/pheno.xlsx"
pheno<- read.xlsx(phenofile, colNames = FALSE, startRow=2)

pheno.pfc <- pheno[which(pheno$X9=="frontal cortex" & pheno$X14 != "braak.stage: Exclude") , c("X1", "X14")]

library(stringr)
pheno.pfc$stage <- str_split_fixed(pheno.pfc$X14, ":", n=2)[,2]

# merge
exp<- mvalues.filtered[, pheno.pfc$X1]

# single cpg p-values

design = model.matrix(~ stage, data = pheno.pfc)

identical(as.character(pheno.pfc$X1), colnames(exp))

#one<- exp[2,]
#coef(summary(lm(as.numeric(one)~as.numeric(stage), data=pheno.pfc)))

lmF=function(methylation){  
  tmp=coef(summary(lm(methylation ~ as.numeric(stage), data=pheno.pfc)))
  return(tmp[2,]) 
}

res <- t(apply(exp,1, lmF))

write.csv (res, file="single-cpg-pvalues.csv")


#result_dms <- cbind (res, p.adjust(res[ ,4], method="fdr") )
#colnames(result_dms)[5] <-"fdr"

