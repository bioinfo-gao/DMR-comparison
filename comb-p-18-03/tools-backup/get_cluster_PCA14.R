#install.packages("FactoMineR")
library(FactoMineR)
library(stringr)

beta_value=readRDS("../betaValue.rds")                       #####-----#### input your own beta value                        
m_value <- log2(beta_value/(1-beta_value))
annot_450K_Aclust=readRDS("../annot-Aclust-format.rds")


mvalue_for_subset_training <- vector("list",50)
mvalue_for_subset_testing <- vector("list",50)

for(i in 1:10){ #### 1
    #i =1
    #cluster_variable_name <- paste("cluster", i, sep = "_")
    #cluster_variable_name_with_5_or_more_CpG <- paste("cluster", i, "with_5_or_more_CpG", sep = "_")
    
    rds_file <- paste0("/media/2T_Disk/Dropbox/ZG_for_LW_Projs/For-Lily/2017_03/complete_dms/dms_pvals_PFC_02_28_", i, ".rds")
    
    tmp <- readRDS(rds_file)
    #assign(cluster_variable_name, tmp)    
    tmp[[1]] <- tmp[[1]][lapply(tmp[[1]],length) > 4] #### 1 means chromsome 1
    #assign(cluster_variable_name_with_5_or_more_CpG, tmp[[1]])  # contain 31 list

    patients_50_selection <- read.csv("50-times-sampling1.csv", row.names=1,header=T)
    
    k <- i+2 
    patients_training <-  patients_50_selection[patients_50_selection[,k] == T,]
    mvalue_for_subset_training[[i]] <-  m_value[colnames(m_value) %in% row.names(patients_training)] #?? changed from patients_training[i]
    patient_stage <- patients_training[2]
    list_of_clusters <- list(data.frame()) # creat a list of matrix
    for( j in 1:length(tmp[[1]]) ) {
        list_of_clusters[[j]] <- mvalue_for_subset_training[[1]][row.names(mvalue_for_subset_training[[1]]) %in% tmp[[1]][[j]],] # chrom 1
        
        t_m_value <- t(list_of_clusters[[j]]) ## transpose of m_value
        patient_mvalue_stage <- merge(t_m_value, patient_stage , by = 0) # by =0  or by = "row.names"
        
        stage_category <- patient_mvalue_stage[,ncol(patient_mvalue_stage)]
        #stage_category <- sapply(stage_category, str_sub(-2, -1))
        stage_category <- str_sub(stage_category, -1, -1)
        
        
        #stage_category <- as.factor(patient_mvalue_stage[,ncol(patient_mvalue_stage)])
        pdf(paste0("PCA/FactMineR_",i,"_cluster_",j,".pdf"),width=15,height=15) #### for the big 45 x 45 matrix
        res.pca = PCA(patient_mvalue_stage[,-c(1,ncol(patient_mvalue_stage))], scale.unit=TRUE, ncp=10, graph=T)#),
        dev.off()
        
        pdf(paste0("PCA/Color_",i,"_cluster_",j,".pdf"),width=15,height=15) #### for the big 45 x 45 matrix
        plot(res.pca$ind$coord[,1], res.pca$ind$coord[,2], col=stage_category) ###
        #text(res.pca$ind$coord[,1], res.pca$ind$coord[,2], rownames(patient_mvalue_stage), pos= 3 )
        text(res.pca$ind$coord[,1], res.pca$ind$coord[,2], stage_category, pos= 3 )
        dev.off()
        
        
        pdf(paste0("PCA/Variance_Cycle_",i,"_cluster_",j,".pdf"),width=5,height=5) #### for the big 45 x 45 matrix
        plot(res.pca$eig[,2], ylab = "Percent of Variance")
        dev.off()
        
        pdf(paste0("PCA/Accumulated_Variance_Cycle_",i,"_cluster_",j,".pdf"),width=5,height=5) #### for the big 45 x 45 matrix
        plot(res.pca$eig[,3], ylab = "Percent of Variance")
        dev.off()
    }

}


#rm(list=ls())

