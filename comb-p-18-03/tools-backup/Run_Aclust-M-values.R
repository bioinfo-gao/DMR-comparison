install.packages("devtools") 
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite() #install automatically 'Biobase' 'IRanges' 'AnnotationDbi' ‘BiocGenerics’ ‘RSQLite’
install_github("tamartsi/Aclust", dependencies=TRUE)
library(Aclust)

setwd("/media/2T-Disk/Dropbox (BBC)/Zhen_Gao/For-Lily")   #####-----#### set your own wd 

beta_value=readRDS("betaValue.rds")                       #####-----#### input your own beta value                        
m_value <- log2(beta_value/(1-beta_value))
annot_450K_Aclust=readRDS("annot-Aclust-format.rds")

chrome_annot_files <- vector("list",24) 
clusters_m_value_list_files <- vector("list",24) 
chromosome_list <- c(1:22, "X", "Y")

#### Apply the Aclust to do clusting , have to do chromosome one by one
for (i in 1:24){
    chromosome <- chromosome_list[[i]]
    chrome_annot_files[[i]] <- annot_450K_Aclust[(annot_450K_Aclust$CHR == chromosome), ]
    chrome_annot_files[[i]] <- subset(chrome_annot_files[[i]], !duplicated(chrome_annot_files[[i]][,4]))       #### Aclust don't allow duplicate in chromsome locus
    clusters_m_value_list_files[[i]] <- assign.to.clusters(m_value, chrome_annot_files[[i]], dist.thresh = 0.5, bp.merge = 1000) # clustring based on m-value
}

#### print out the clusters
for (i in 1:24){
    count <- 0
    chromosome <- chromosome_list[[i]]
    print(paste0("For ", chromosome, ", there are following DMR" ))
    for (item in clusters_m_value_list_files[[i]]){
        if (length(item)>4) {  # choose the cluster contain 5 or more CpG site
            count <- count + 1 
            df <- chrome_annot_files[[i]][chrome_annot_files[[i]]$IlmnID %in% item]  # get the annotation in that cluster
            start_pos <- head(df, 1)[,4]
            end_pos   <- tail(df, 1)[,4]
            region_info <-paste0("chromosome:", chromosome, "-", start_pos[[1]], "-", end_pos[[1]])
            CpG_and_locus <-  df[,c(1,4)]
            order <- rep(count, length(item))
            region <- rep(region_info, length(item))
            CpG_in_this_cluster <- rep(length(item), length(item))
            output <- cbind(order, CpG_and_locus, region, CpG_in_this_cluster)
            print(output) 
        }
    }
}

##### initiate the matrix to store all mean bata-values for 110 patients
mean_rows_m_value_in_cluster <- matrix( , nrow = 0, ncol = ncol(beta_value)) #### m-value matrix, clumne equals the patient number
mean_rows_m_value_in_cluster_for_each_person <- vector(mode="numeric", length=ncol(beta_value))  # the m-value vector for each person, each vector for a specific CpG cluster

each_cluster_m_value_vector_list <- list()           # A list of vectors,  the number of vectors is the cluster numbers, each vector contain the mean m-value for that cluster of all patients
cluster_names <- vector(mode="character", length=0)  # each cluster has different name

all_names <- vector(mode="character", length=0)  # the names for a clusters, start from chromsome 1 to chromsome Y, eg. Chromome-X-CpG-cluster-235
all_count <- 0
for (i in 1:24){
    count <- 0
    for (item in clusters_m_value_list_files[[i]]){
        if (length(item)>4) {  # calculate the mean value if the cluster contains 5 or more CpG sites
            count <- count + 1
            all_count <- all_count + 1
            for (p in 1:ncol(m_value)) {                # loop from the first patient to last patient
                Patient_CpG_value <-  m_value[p]            # get the m value belong to certain person 
                mean_rows_m_value_in_cluster_for_each_person[p] <- mean(Patient_CpG_value[row.names(Patient_CpG_value) %in% item,])  # get the mean m value in that cluster
            }
            chromosome <- chromosome_list[[i]]
            cluster_names[count] <- paste0(chromosome, "-CpG-cluster-", toString(count)) ##### e.g. Chromome-X-CpG-cluster-235, count from 0 for each chrom
            all_names[all_count] <- paste0(chromosome, "-CpG-cluster-", toString(count)) ##### e.g. Chromome-X-CpG-cluster-235 but count from chrom 1 to chrom Y, used outside the for loop
            values <- mean_rows_m_value_in_cluster_for_each_person #####
            each_cluster_m_value_vector_list[[cluster_names[count]]] <- values
            mean_rows_m_value_in_cluster <- rbind(mean_rows_m_value_in_cluster, each_cluster_m_value_vector_list[[cluster_names[count]]]) 
        }
    }
}

rownames(mean_rows_m_value_in_cluster) <- all_names # give all clusters a specific name
write.csv(mean_rows_m_value_in_cluster, "mean_rows_beta_value_in_cluster.csv")
#mean_in_cluster_and_stage <- rbind(mean_rows_m_value_in_cluster, stage) 
#write.csv(mean_in_cluster_and_stage, "mean_m_in_cluster_and_stage.csv")
