library(glmnet) #
library(caret)
library(dplyr)
library(pROC)

#########   1. prepare some data structure and read beta value 
Accuracy <- vector()    # for split value from train sample  #
Sensitivity <- vector()                                      #
Specificity <- vector()                                      #
AUC <- vector()                                              #
Accuracy1 <- vector()   # for split value from train sample  #
Sensitivity1 <- vector()                                     #
Specificity1 <- vector()                                     #
AUC1 <- vector()                                             #

#########   1. prepare some data structure and read beta value 
setwd("~/2T_Disk/Dropbox/Zhen_Gao/LW_Projects/For-Lily/2017-09/2017-09-28-PC1-selection-repeated_PC1-prediction")
data.dir <- "/media/2T_Disk/Dropbox/Zhen_Gao/LW_Projects/For-Lily/2017-08/Annot/"
beta_value=readRDS(paste0(data.dir, "betaValue.rds"))                           
annot_450K=readRDS(paste0(data.dir, "annot_450K_chr_coord.rds")) # no sex chromsome, non-specific     
#m_value_all_patients <- get_m_value_wo_filter_of_low_variance(beta_value, annot_450K) ########### 5. filter according to methylation level 
m_value_all_patients <- get_m_value_with_filter_of_low_variance(beta_value, annot_450K) ########### 5. filter according to methylation level 

patients <- read.csv("../sorted-AD-stage1234-stage56.csv", row.names=1)
nfold.cv <- 5

for (n in 1:10) {
    #n=1
    seed = 122 + n 
    set.seed(seed)
    flds <- createFolds(patients$Stage, k = nfold.cv, list = TRUE, returnTrain = F)    
    
    for( i in 1:nfold.cv ){ #### 1
        #i = 1  #########   
        patients.train = patients[-(flds[[i]]),]
        patients.test  = patients[flds[[i]],]
        #train_PC1_component <- list()
        #test_PC1_component <- list()
        training_PC1_matrix <- list()
        testing_PC1_matrix <- list()
        #rds_file <- paste0("./cor_cluster/highly_correlated_list_", i, ".rds")#HRER
        rds_file <- paste0("./cor_cluster/highly_correlated_list_low_variance_filtered", i, ".rds") #HRER match the filter or not above
        total_cluster <- readRDS(rds_file) #oooooo
        
        m_value_train <- m_value_all_patients[,match(rownames(patients.train), colnames(m_value_all_patients))]  
        m_value_train_ordered <- m_value_train[match(rownames(annot_450K), rownames(m_value_train)),] 
        m_value_test <- m_value_all_patients[,match(rownames(patients.test), colnames(m_value_all_patients))]  
        m_value_test_ordered <- m_value_test[match(rownames(annot_450K), rownames(m_value_test)),] 
        
        # res.list.long <- total_cluster [lapply(total_cluster, length)>=3 ]
        # cluster.table <- cpglist.to.table (cpgs.list=res.list.long, methylvalue = beta_value)
        
        # cluster_mean <- aggregate(x = cluster.table, by = list(unique.values = cluster.table$cluster), FUN = mean)
        # cluster_mean_value <- cluster_mean[,-c(1,2)]
        # cluster_mean_value_ordered <-  cluster_mean_value[,match(rownames(patients.train), colnames(cluster_mean_value))] 
        # cluster_mean_value_ordered <-  as.matrix( cluster_mean_value_ordered)  
        
        train_patient_stage <-  patients.train$Stage 
        
        #for( k in 1:nrow(cluster_mean_value_ordered)) {          
        for( k in 1:length(total_cluster)) {          
            #k=1
            training_cluster <- m_value_train_ordered[row.names( m_value_train_ordered ) %in% total_cluster[[k]],]
            testing_cluster  <- m_value_test_ordered[row.names( m_value_test_ordered ) %in% total_cluster[[k]],]
            
            training_cluster_pca = prcomp(t(training_cluster), center=T, scale. = T)
            train.data <- as.data.frame(training_cluster_pca$x) # PCs, see genomicsclass.github.io/book/pages/pca_svd.html
            
            # train_cluster_mean <- cluster_mean_value_ordered[c(k),]
            ttest_result  <- t.test( train.data$PC1 ~ train_patient_stage) 
            #print(ttest_result$p.value)
            print(k)   
            
            if(ttest_result$p.value < 0.05) {
                
                test.data <- predict(training_cluster_pca, newdata = t(testing_cluster))
                test.data  <- as.data.frame(test.data)
                
                print("***********")
                print(ttest_result$p.value)
                
                training_PC1_matrix <- cbind(training_PC1_matrix, train.data$PC1)
                testing_PC1_matrix <- cbind(testing_PC1_matrix, test.data$PC1)
            }
            
        }
        
        #train_PC1_component[[k]] <- data.frame(training_cluster_pca$x)[,1, drop = FALSE] # PC1 see genomicsclass.github.io/book/pages/pca_svd.html
        #test_PC1_component[[k]] <- predict(training_cluster_pca, newdata = t(testing_cluster))[,1,drop = FALSE] #transform test into PCA then get PC1
        #}
        
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(i)
        
        split.value <- get_split_value_from_train_data(patients.train)
        Stage <- patients.train$Stage
        
        training_PC1_matrix <-  apply(training_PC1_matrix,2, unlist) 
        testing_PC1_matrix  <-  apply(testing_PC1_matrix,2, unlist) 
        
        train_patient_stage_and_PC1s = as.data.frame(cbind(Stage,as.data.frame(training_PC1_matrix)))
        test_patient_PC1s = as.data.frame(testing_PC1_matrix)
        
        colnames(train_patient_stage_and_PC1s) = append("Stage", paste0("PC1_", 1:ncol(test_patient_PC1s)))
        rownames(train_patient_stage_and_PC1s) = rownames(patients.train)
        
        colnames( test_patient_PC1s) = paste0("PC1_", 1:ncol(test_patient_PC1s))
        rownames( test_patient_PC1s) = rownames(patients.test)
        
        levels(train_patient_stage_and_PC1s$Stage) <- make.names(levels(factor(train_patient_stage_and_PC1s$Stage))) #levels are needed for train()fun
        
        # https://stats.stackexchange.com/questions/250187/variable-selection-with-elastic-net-with-r-glmnet-carret
        # https://rstudio-pubs-static.s3.amazonaws.com/43302_2d242dbea93b46c98ed60f6ac8c62edf.html
        lambda.grid <- seq(0, 100)
        alpha.grid <- seq(0, 1.0, length = 11)
        srchGrd = expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)
        fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, classProbs = TRUE, summaryFunction = twoClassSummary) #repeated CV
        
        eNetModel <- train(Stage ~ ., 
                           #data=train_patient_stage_and_mvalue, 
                           data=train_patient_stage_and_PC1s, 
                           method = "glmnet", 
                           metric="ROC", 
                           trControl = fitControl, 
                           family="binomial", 
                           tuneGrid = srchGrd) #,tuneLength=10
        
        pred.eNetModel <- predict(eNetModel, test_patient_PC1s, type="prob")[,"Late.stage"]#type must be either "raw" or "prob"
        #pred.eNetModel <- predict(eNetModel, train_patient_stage_and_mvalue[,-(1)], type="prob")[,"Late.stage"]#type must be either "raw" or "prob"
        
        ROC <- pROC::roc(patients.test$Stage, pred.eNetModel)
        AUC[i] <- AUC1[i] <- pROC::auc(ROC)  
        
        pred.value1 <- pred.value <- as.numeric(pred.eNetModel)   
        
        # very veirld mistake cannot determine 7.88239176324881e-05 < 0.5
        #if (pred.value < split.value)  {pred.value <- "Before-late"}
        #else { pred.value <- "Late-stage" }  # this one is much better than split 0.5
        
        pred.value[ as.numeric(pred.value1) >= split.value ] <- "Late-stage"
        pred.value[ as.numeric(pred.value1) < split.value ]  <- "Before-late" 
        
        my.confusionMatrix<- confusionMatrix(pred.value, patients.test$Stage)
        
        Sensitivity[i] <- my.confusionMatrix$byClass["Sensitivity"]
        Specificity[i] <- my.confusionMatrix$byClass["Specificity"]
        Accuracy[i]    <- my.confusionMatrix$overall["Accuracy"]
        
        ####
        pred.value1[ as.numeric(pred.value1) >= 0.5 ] <- "Late-stage"
        pred.value1[ as.numeric(pred.value1) < 0.5 ]  <- "Before-late" 
        
        my.confusionMatrix1<- confusionMatrix(pred.value1, patients.test$Stage)
        
        Sensitivity1[i] <- my.confusionMatrix1$byClass["Sensitivity"]
        Specificity1[i] <- my.confusionMatrix1$byClass["Specificity"]
        Accuracy1[i]    <- my.confusionMatrix1$overall["Accuracy"]
    }   
    
    general_performance_results <- cbind(Sensitivity, Specificity, Accuracy, AUC)
    prediction_performance_file <- paste0("./PC1/filter/" ,n, "/","correlation_cluster_mean_cv_09_28_train_split", ".csv")
    write.csv(general_performance_results, prediction_performance_file)
    
    general_performance_results1 <- cbind(Sensitivity1, Specificity1, Accuracy1, AUC1)
    prediction_performance_file1 <- paste0("./PC1/filter/", n, "/", "correlation_cluster_mean_cv_09_28_point_5_split", ".csv")
    write.csv(general_performance_results1, prediction_performance_file1)
}