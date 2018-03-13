source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("caTools")
library(caTools)

data.dir <- "/home/gao/R/Lily_PFC"
setwd(data.dir)
patients <- read.csv(paste0(data.dir,"/sorted1.csv"), row.names=1)
smp_size <- floor(0.5 * nrow(patients))  ## 50% of the sample size
set.seed(123)                            ## set the seed for partition reproductivity

for (val in 1:50) {
  train_rows = sample.split(patients$Stage, SplitRatio=0.5)
  patients <- cbind(patients,train_rows)
}

cols <- sapply(patients, is.logical)
patients[,cols] <- lapply(patients[,cols], as.numeric)

write.csv(patients, file = "50-times-sampling1.csv")
