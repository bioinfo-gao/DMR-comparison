#!/usr/bin/Rscript
# authored by ZG
# summary: identify Differentially methylated clusters with Aclust in 450K data

library(devtools)
library(Aclust)

# args = commandArgs(trailingOnly=T)
# 
# i = args[1]
# #i = argsL$i1
# i <- as.integer(i) 

# dir <- "/nethome/zxg161/For-Lily"
# setwd(dir)

#### 1.1 readin the beta value and change to M-value 
beta_value=readRDS("betaValue.rds")                       #####-----#### input your own beta value                        
m_value <- log2(beta_value/(1-beta_value))
annot_450K_Aclust=readRDS("annot-Aclust-format.rds")

## get subset of m_value and subset of PFC
patients_50_selection <- read.csv(paste0(dir,"/50-times-sampling1.csv"), row.names=1,header=T)
patients_50_selection <- patients_50_selection[-c(1:2)]

mvalue_for_subset_training <- vector("list",50) 
mvalue_for_subset_testing <- vector("list",50) 


Run_Aclus <- function(mval, distance_method = "complete") {
    chrome_annot_files <- vector("list",24) 
    clusters_m_value_list_files <- vector("list",24) 
    chromosome_list <- c(1:22, "X", "Y")
    for (j in 1){ #### change from 1:24
        chromosome <- chromosome_list[[j]]
        chrome_annot_files[[j]] <- annot_450K_Aclust[(annot_450K_Aclust$CHR == chromosome), ]
        chrome_annot_files[[j]] <- subset(chrome_annot_files[[j]], !duplicated(chrome_annot_files[[j]][,4]))       #### Aclust don't allow duplicate in chromsome locus
        clusters_m_value_list_files[[j]] <- assign.to.clusters(m_value, chrome_annot_files[[j]], dist.type = "pearson", method = distance_method, dist.thresh = 0.5, bp.merge = 1000) # clustring based on m-value
    }
    return (clusters_m_value_list_files)
}

patients_training <-  patients_50_selection[patients_50_selection[,i] == T,]
mvalue_for_subset_training[[i]] <-  m_value[colnames(m_value) %in% row.names(patients_training)] # changed from patients_training[i] not very sure
    
result_dms_training <- Run_Aclus(mvalue_for_subset_training[[i]])
training_outfile <- paste0("dms_pvals_PFC_02_28_", toString(i), ".rds")
saveRDS(result_dms_training, training_outfile)
