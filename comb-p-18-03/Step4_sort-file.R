#dir<-"/media/2T-Disk/DMS/p-values/TSV"
#out.dir <- "/media/2T-Disk/DMS/p-values/For-comb-p/"
  
df_train  <- read.table("pvals_PFC_TSV.bed", sep="\t", header=TRUE)
ordered_traing <- df_train[order(df_train[,1],df_train[,2], decreasing=F),]
train_bed_file <- paste0("Ordered_Training.bed")
write.table (ordered_traing, train_bed_file, sep="\t",row.names=FALSE)
    
