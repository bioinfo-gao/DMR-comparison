#install.packages("dplyr")
library("dplyr")
#dir<-"~/2T_Disk/Dropbox/Zhen_Gao/LW_Projects/comb-p-2017-09-rerun"
#setwd(dir)

#### 1.1 read in chromsome and locus info
#locus_file <- read.csv(paste0(dir, "/HumanMethylation_simple1.csv"), header=T)
locus_file <- read.csv("HumanMethylation_simple1.csv", header=T)
vars <- c("chrom", "start", "end", "P.Value")

#train_file_name <- "dms_pvals_PFC_6stage.csv"
train_file_name <- "dms_pvals_PFC_Before_Late_vs_Late.csv"
pvalue_file <- read.csv(train_file_name, header=T)
train_locus_p_value_file <- inner_join(locus_file, pvalue_file, by.locus_file = "X", by.pvalue_file = "X")

train_locus_p_value_file$end <- train_locus_p_value_file[["MAPINFO"]]+25 ##########
train_locus_p_value_file$start <- train_locus_p_value_file[["MAPINFO"]]-25 ###########
train_locus_p_value_file$chrom <- train_locus_p_value_file[["CHR"]]
train_locus_p_value_file$P.Value <- train_locus_p_value_file[["fdr"]]

training_outfile <- train_locus_p_value_file[vars]
train_bed_file <- "pvals_PFC_training.bed" 
write.csv (training_outfile, train_bed_file)
