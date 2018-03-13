#setwd("C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/RealDataAnalysis")
setwd("/media/2T_Disk/Dropbox/Zhen_Gao/Saurav-DMR")
beta.value.all <- readRDS("beta.value.all.realdata.rds")
dim(beta.value.all)#probe=333757, sample=62

mvalues <- log2(beta.value.all/ (1- beta.value.all))

groups <- rep(c("tumore", "normal"), each = 31)

pheno <- data.frame( cbind (colnames(mvalues), groups))

#mval <- mvalues[1 ,]

lmF <- function(mval) {
  tmp = coef(summary(lm ( mval~ groups, data=pheno)))
  return (tmp[2,1:4])
}

allcpg <- data.frame(t(apply(mvalues,1,lmF)))

allcpg$ILMNID <- row.names(allcpg)

saveRDS (allcpg, "singleCpG.pvals.RDS")

############# add annotations

cpg.location <- readRDS ("C:/Users/lxw391/Box Sync/METHODS-SHARED/METHOD-METHYL-GENE-TEST/AD/DATA/cpg.locations.RDS" )

result.location <- merge (allcpg, cpg.location, by =c("ILMNID", "id"))

result.location$chrom <- paste0("chr", result.location$CHR)
result.location$start <- as.numeric(as.character(result.location$MAPINFO)) - 1
result.location$end   <- result.location$MAPINFO

result.location$PVALS <- result.location$`Pr...t..`

result.location <- result.location[order(result.location$CHR, as.numeric(as.character(result.location$MAPINFO))) ,]

result.final   <- result.location[, c("chrom", "start", "end", "PVALS")]

write.csv(result.final,
          file=paste0("linear_model_result_real_data.csv"), 
          row.names = FALSE)

