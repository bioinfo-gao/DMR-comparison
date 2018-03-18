setwd("D:/Saurav SXM1600 backup/R Script files") #-------- need access to this file------------- #

#setwd("C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017")

###debug(utils:::unpackPkgZip)#######use this before using "install.packages" command if "unable to move temporary..." error comes
###install.packages("xml2")

#source("https://bioconductor.org/biocLite.R")
#biocLite("RnBeads")
#library("RnBeads")

#source("https://bioconductor.org/biocLite.R")
#biocLite("BayesPeak")
#biocLite("genomation")
library("genomation")
library(IRanges)
library(GenomicRanges)

# tt<-scan(filenames[5],sep=',',what=list('character','character','character','numeric','numeric','numeric','numeric'))
# treatment <- read.bed(tt)
#treatment <- read.bed(filenames[5])##??

# myfile<-system.file(paste0("D:/Saurav SXM1600 backup/R Script files/Comb-P-Folders/linear_model_result_for_mu",miuv[8],"_repetition_",rp[5]),"Train.regions-t.bed",package="genomation")
# refseq1<-system.file(myfile,package="genomation")
# te<-readBed(refseq1, track.line = FALSE, remove.unusual = FALSE,
#         zero.based = TRUE)

# dir <- system.file("extdata", package="BayesPeak")
# file <- file.path(dir, "H3K4me3reduced.bed")
# treatment <- read.bed(file)
# 
# my.file=system.file("extdata","chr21.refseq.hg19.bed",package="genomation")
# refseq = readBed(my.file,track.line=FALSE,remove.unusual=FALSE)
# head(refseq)
# 
# treatment
# data.dir <- "C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017/data/"
# data.dir1 <- "D:/Saurav SXM1600 backup/R Script files/Comb-P-Folders/"  # ----- ?????----- #
# result.dir <- "C:/Users/lxw391/Dropbox (BBSR)/Saurav (Chu) second/DMR-comparions/LWtest/Test_11-10-2017/results_bumphunter_1000permutation/"
## not used data.dir <- "/media/2T_Disk/Dropbox (BBSR)/Saurav_(Chu)_second/DMR-comparions/LWtest/Test_11-10-2017/data/"
data.dir1 <- "D:/Saurav SXM1600 backup/R Script files/Comb-P-Folders/"  # ----- ?????----- #
result.dir <- "/media/2T_Disk/Dropbox (BBSR)/Saurav_(Chu)_second/DMR-comparions/LWtest/Test_11-10-2017/results_bumphunter_1000permutation/"


############### 3. annotation files -------------------------------------------------
load("/media/2T_Disk/Dropbox/Zhen_Gao/MethylationFunctions/fullannotInd.rda")   # data downloaded from https://rforge.net/IMA/ , under "4 Annotation file"
annot <- as.data.frame(fullannot)
vars <- c("ILMNID", "CHR", "MAPINFO")
cpg.location <- annot[vars]


######Combp result extraction######
seed.values <- c(100, 210, 330, 450, 680) ##select random clusters from A-clusters 

miu_values <- c(0, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40) # simulation scenarios

cpg_number_cutoffval <- 4 ##newly added #select those clusters with >4 cpgs


# 5. run Comb-P method
source("3_runCombp_myLW2LW2.R") ##corrected the name

result.all = as.data.frame(matrix(ncol=10, nrow=length(seed.values)*length(miu_values)))##newly added

for (m in 1:length(miu_values)) ##newly added
{   #m=1
    print(paste0("miu_values=",miu_values[m]))
    
    for (r in 1:length(seed.values)) ##newly added
    {
        #m=8
        #r=2
        
        # 2. run comb-p and process results
        mval_each <- miu_values[m]##newly added
        temp <- readRDS(paste0("simData_miu",miu_values[m],"_rep",r,".rds") )
        
        # merge cpgs of beta value data with whole-genome wide annotation file to obtain whole-genome wide cpg's beta value
        temp$cpg <- rownames(temp)#newly added
        
        temp.location <- merge(temp, cpg.location, by.x="cpg", by.y="ILMNID")
        rownames(temp.location) <- temp.location$cpg
        temp.location$CHR <- paste0("chr", temp.location$CHR)
        temp.location$MAPINFO <- as.numeric(as.character(temp.location$MAPINFO))
        
        #temp.location <- temp.location[order(temp.location$CHR, temp.location$MAPINFO) ,]
        
        location.index <- c(grep("CHR", colnames(temp.location)), grep("MAPINFO", colnames(temp.location)), grep("cpg", colnames(temp.location)) )
        
        beta.matrix <- as.matrix (temp.location[, -location.index])  #take only the annotated beta values
        
        ##type <- rep(c("cancer","normal"),c(exsampleno,ctsampleno))
        ##design1 <- model.matrix(~type)  ##as "~type" is only used, so coef=2 (2 columns) here
        
        #data.dir1 <- "D:/Saurav SXM1600 backup/R Script files/Comb-P-Folders/"
        
        #filenames <- list.files(path =paste0(data.dir1,"linear_model_result_for_mu",mval_each,"_repetition_",r),  
        #                         pattern=".bed", full.names=TRUE)
        
        #combp.result <- read.table(filenames[5],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
        
        filenm=paste0(data.dir1,"linear_model_result_for_mu",mval_each,"_repetition_",r,"/Train.regions-t.bed")
        combp.result <- tryCatch({
            if(file.size(filenm)>0) {
                read.table(filenm, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
            }
        }, error = function(e1){return(NULL)})
        
        # extract results
        if(is.null(combp.result)==FALSE) ####if any predicted cluster is identified from Combp
        {  
            colnames(combp.result) <- c("#chrom","start","end","min_p","n_probes","z_p","z_sidak_p")
            
            #elapsedtime <-  proc.time() - ptm ##newly added the line
            
            result <- combp.result[which(combp.result$z_p < 0.05 & combp.result$n_probes > cpg_number_cutoffval), ]
        }
        
        if(is.null(combp.result)==TRUE) ####if no predicted cluster is identified from Combp ##newly added
        {
            result <- NULL
            
            #elapsedtime <-  proc.time() - ptm ##newly added the line
        }
        
        
        temp.ranges <- result
        #Elapsedtime <- result[[2]] ##??
        Elapsedtime <- "Not found"
        
        # process results    
        
        simAclust.file <- read.csv (paste0("simAclust_for",miu_values[m], "_rep" , r, ".csv"))
        
        one <-  processCombpResults (combp.ranges=temp.ranges, aclustfile=simAclust.file, mval_each, r, cpg_number_cutoffval)##updated
        
        one$elapsedtime <- Elapsedtime ##newly added
        one$mu <- miu_values[m]
        one$rep <- r
        one$method <- "Combp" ##newly added
        
        names(result.all) <- names(one) ##newly added
        result.all[((m-1)*length(seed.values)+r),] <- one ##reformated
        write.table(result.all,file=paste("Combp_result.all_for_",length(miu_values),"_number_of_mu_values_having_repitition_of_",length(seed.values),"_times.csv",sep=""),sep=",",row.names=F,col.names=T)##newly added
        
    }
    
}














