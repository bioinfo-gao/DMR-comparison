# main file
setwd("C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017/files-for-combP")

data.dir <- "C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017/data/"



################# 1. Beta values for all cpgs in the genome ---------------------------------
# note this file is the same for every simulation scenario

beta.value.all <- readRDS(paste0(data.dir, "beta.value.all.rds"))

exsampleno<-7 ##the number of samples in the First group
ctsampleno<-7 ##the number of samples in Second group

colnames(beta.value.all)[1:(exsampleno)]<-paste(colnames(beta.value.all)[1:(exsampleno)],'-','Tumour',sep = "")
colnames(beta.value.all)[(1+exsampleno):ncol(beta.value.all)] <-paste(colnames(beta.value.all)[(1+exsampleno): ncol(beta.value.all)],'-','Normal',sep = "")

################# 2. files for all scenarios ---------------------------------------------- 
res_aclcpg=read.csv(paste0(data.dir, "A-clust-results.csv"),header=T,sep=",") ## with cpgnames (rows) and samples (columns)###result of a-clustering

colm<-ncol(res_aclcpg)

beta_start_columnid<-9 ##newly added
  
colnames(res_aclcpg)[beta_start_columnid:(beta_start_columnid+exsampleno-1)]<-paste(colnames(res_aclcpg)[beta_start_columnid:(beta_start_columnid+exsampleno-1)],'-','Tumour',sep = "")
colnames(res_aclcpg)[(beta_start_columnid+exsampleno):colm]<-paste(colnames(res_aclcpg)[(beta_start_columnid+exsampleno):colm],'-','Normal',sep = "")
rownames(res_aclcpg)<-res_aclcpg$cpg

totalclsno_aclst<-max(res_aclcpg$Clusternumber) ###find the total number of clusters initially from Aclust results

seed.values <- c(100,210,330,450,680) ##select random clusters from A-clusters 

m.values <- c(0,  0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40) # simulation scenarios


# seed.values <- c(100,210)
# 
# m.values <- c(0.3, 0.4) # simulation scenarios


cpg_number_cutoffval<-4 ##newly added #select those clusters with >=5 cpgs


############### 3. annotation files -------------------------------------------------
load (paste0(data.dir, "fullannotInd.rda"))   # data downloaded from https://rforge.net/IMA/ , under "4 Annotation file"
annot <- as.data.frame(fullannot)
vars <- c("ILMNID", "CHR", "MAPINFO")
cpg.location <- annot[vars]
cpg.location$ILMNID <- as.character(cpg.location$ILMNID)

# 2. run linear model

runLinearModel<-function(beta.values, cpg_number_cutoff, mval_each, r) ##newly added last three arguments
{  
  ptm <- proc.time() ##newly added the line
  require(DMRcate)
  myBetas<-as.matrix(beta.values)
  myMs <- logit2(myBetas)  ##log base 2 transformation
  
  type <- factor(sub(".*-", "", colnames(myMs)))
  design1 <- model.matrix(~type)  ##as "~type" is only used, so coef=2 (2 columns) here
  
  # one<- myMs[2,]
  # coef(summary(lm(as.numeric(one) ~ type))) [2 ,]
  
  #res <- as.data.frame(matrix(ncol=4, nrow=0))

  #   for (i in 1:nrow(myMs)){
  #   print (i)
  #   tmp=coef(summary(lm(myMs[i ,] ~ type)))
  #   colnames(res) <- colnames(tmp)
  #   res <- rbind (res, tmp[2,])
  # }
 
  lmF <- function(methylation){
     tmp=coef(summary(lm(methylation ~ type)))
     return(tmp[2,])
   }

  res <- as.data.frame( t(apply(myMs,1, lmF)) )
  res$id <- rownames(res)
  
  return(res)
}



for (m in 8:8) 
{
  
  for (r in 1:length(seed.values)) 
  {
    #m=8
    #r=1

     #2. run DMRcate and process results
     mval_each<-m.values[m]##newly added
     temp <- readRDS (paste0("C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017/simulation-datasets/simData_miu",m.values[m],"_rep",r,".rds") )
    
     result <- runLinearModel(beta.values=temp,
                              cpg_number_cutoff=cpg_number_cutoffval,
                              mval_each, r)
    
     result.location <- merge (result, cpg.location, by.x="id", by.y="ILMNID")
     
     result.location$chrom <- paste0("chr", result.location$CHR)
     result.location$start <- as.numeric(as.character(result.location$MAPINFO)) - 1
     result.location$end   <- result.location$MAPINFO
     
     result.location$PVALS <- result.location$`Pr(>|t|)`
     
     result.location <- result.location[order(result.location$CHR, result.location$MAPINFO) ,]
    
     result.final   <- result.location[, c("chrom", "start", "end", "PVALS")]
     
     write.csv(result.final,
              file=paste0("linear_model_result_for_mu",mval_each,"_repetition_", r,".csv"), 
              row.names = FALSE)
    
  }
}

