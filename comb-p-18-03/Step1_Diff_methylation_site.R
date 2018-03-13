# authored by LW and ZG
# summary: identify Differentially methylated sites in 450K data
#instalmF <- function(mval, subPFC) {
#tmp = coef(summary(lm(mval~ stage, data=subPFC)))
#return (tmp[2,1:4])
#ll.packages("MHTdiscrete")
library(MHTdiscrete)

setwd("/media/2T_Disk/Dropbox/Zhen_Gao/Saurav/comb-p-18-03")

#data.dir <- "/home/gao/H-Driver-Link/Zhen-Gao/ExampleCode/DMS/Source-Data"
data.dir <- "/media/2T_Disk/Dropbox/HD_Zhen_Gao/ExampleCode/DMS/Source-Data"
#dir<-"~/2T_Disk/Dropbox/Zhen_Gao/LW_Projects/comb-p-2017-09-rerun"
#setwd(dir)

#### 1.1 read in phenotype info
pheno <- read.csv(paste0(data.dir, "/pheno.csv"), header=FALSE)
PFC <- pheno[which(pheno$V8=="frontal cortex" & pheno$V13 !="braak.stage: Exclude") ,] # & pheno$V12 != "ad.disease.status: Exclude"
#retained only row 2 as headers
betas <- read.csv(paste0(data.dir,"/GSE59685_betas-PFC2-passed-QC.csv")) 
beta <- betas[, as.character(PFC$V2)]
rownames(beta) <- betas$cpg # or CpG , need take a look to check

mvalue <- log2(beta/(1-beta))

### 2. extract age.Brain, batch(methylation_plate), stage info
PFC$age.Brain <- as.numeric(substring(PFC$V16, 11,14))
PFC$stage <- as.numeric(substring(PFC$V13, 13,14))
# optional but I believe important
PFC$stage <- ifelse(PFC$stage <=4, "Before_Late", "Late_Stage") # shrink 6 stages to 2 stages

PFC$batch <- as.character(substring(PFC$V11, 10,19))
rownames(PFC)<- PFC$V2

vars <- c("age.Brain", "stage", "batch")
PFC <- PFC[vars]

# ### 3. fit linear model 
identical(colnames(mvalue), rownames(PFC) ) #check if sample ids are the same
#identical(colnames(mvalue), rownames(cycle_info)) #check if sample ids are the same


lmF <- function(mval, subPFC) {
    tmp = coef(summary(lm(mval~ stage, data=subPFC)))
    return (tmp[2,1:4])
}

allcpg_training <- t(apply(mvalue,1,lmF,PFC)) 
allcpg_Tra <- allcpg_training 
result_dms <- cbind(allcpg_training, Sidak.p.adjust(allcpg_Tra[ ,4]) ) # Sidak.p.adjust()
colnames(result_dms)[5]<-"fdr"

training_outfile <- paste("dms_pvals_PFC_", "Before_Late_vs_Late.csv", sep="")
write.csv (result_dms, training_outfile)
