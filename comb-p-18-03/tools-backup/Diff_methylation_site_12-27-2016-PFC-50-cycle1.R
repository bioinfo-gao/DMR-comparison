# authored by LW and ZG
# summary: identify Differentially methylated sites in 450K data

data.dir <- "/home/gao/H-Driver-Link/Zhen-Gao/ExampleCode/DMS/Source-Data"
dir<-"/home/gao/2T-Disk/DMS"
setwd(dir)

#### 1.1 read in phenotype info
pheno <- read.csv(paste0(data.dir, "/pheno.csv"), header=FALSE)
PFC <- pheno[which(pheno$V8=="frontal cortex" & pheno$V13 !="braak.stage: Exclude") ,] # & pheno$V12 != "ad.disease.status: Exclude"
#retained only row 2 as headers
betas <- read.csv(paste0(data.dir,"/GSE59685_betas-PFC2-passed-QC.csv")) 
beta <- betas [, as.character(PFC$V2)]
rownames(beta) <- betas$cpg # or CpG , need take a look to check

mvalue <- log2(beta/(1-beta))

#### 1.2 read 50 cycle data
cycle_info <- read.csv(paste0(data.dir, "/50-times-sampling2.csv"), header=T)
cycle_info <- cycle_info[-c(2:4)]
rownames(cycle_info) <- cycle_info$X # or CpG , need take a look to check
cycle_info <- cycle_info[-1] 

### 2. extract age.Brain, batch(methylation_plate), stage info
PFC$age.Brain <- as.numeric(substring(PFC$V16, 11,14))
PFC$stage <- as.numeric(substring(PFC$V13, 13,14))
PFC$batch <- as.character(substring(PFC$V11, 10,19))
rownames(PFC)<- PFC$V2

vars <- c("age.Brain", "stage", "batch")
PFC <- PFC[vars]

# ### 3. fit linear model 
identical(colnames(mvalue), rownames(PFC) ) #check if sample ids are the same
identical(colnames(mvalue), rownames(cycle_info)) #check if sample ids are the same

## get subset of mvalue and subset of PFC
patients_50_selection <- read.csv(paste0(data.dir,"/50-times-sampling1.csv"), row.names=1,header=T)
patients_50_selection <- patients_50_selection[-c(1:3)]

mvalue_for_subset_training <- vector("list",50) 
mvalue_for_subset_testing <- vector("list",50) 
subPFC_training <- vector("list",50) 
subPFC_testing <- vector("list",50) 
allcpg_training <- vector("list",50) 
allcpg_testing <- vector("list",50) 
result_dms_training <- vector("list",50) 
result_dms_testing <- vector("list",50) 

lmF <- function(mval, subPFC) {
    tmp = coef(summary(lm(mval~ stage, data=subPFC)))
    return (tmp[2,1:4])
}

for (i in 1:50) {
    patients_training <-  patients_50_selection[patients_50_selection[,i] == T,]
    patients_testing <-  patients_50_selection[patients_50_selection[,i] == F,]
    mvalue_for_subset_training[[i]] <-  mvalue[colnames(mvalue) %in% row.names(patients_training[i])]
    mvalue_for_subset_testing[[i]] <-  mvalue[colnames(mvalue) %in% row.names(patients_testing[i])]
    subPFC_training[[i]] <- PFC[row.names(PFC) %in% row.names(patients_training[i]),]
    subPFC_testing[[i]] <- PFC[row.names(PFC) %in% row.names(patients_testing[i]),]
    
    allcpg_training[[i]] <- t(apply(mvalue_for_subset_training[[i]],1,lmF,subPFC_training[[i]])) 
    allcpg_Tra <-allcpg_training[[i]] 
    result_dms_training[[i]]<-cbind (allcpg_training[[i]], p.adjust(allcpg_Tra[ ,4], method="fdr") )
    colnames(result_dms_training[[i]])[5]<-"fdr"
    training_outfile <- paste("dms_pvals_PFC_", toString(i), "_training.csv", sep="")
    write.csv (result_dms_training[[i]], training_outfile)
    
    allcpg_testing[[i]] <- t(apply(mvalue_for_subset_testing[[i]],1,lmF,subPFC_testing[[i]]))
    allcpg_Tst <-allcpg_testing[[i]] 
    result_dms_testing[[i]]<-cbind (allcpg_testing[[i]], p.adjust(allcpg_Tst[ ,4], method="fdr") )
    colnames(result_dms_testing[[i]])[5]<-"fdr"
    testing_outfile <- paste("dms_pvals_PFC_", toString(i), "_testing.csv", sep="")
    write.csv (result_dms_testing[[i]], testing_outfile)
}




