#setwd("~/2T_Disk/Dropbox/Tool-Utility/tools/python-tools/comb-p/comb-p-2017-09-rerun-single-CpG")
setwd("/media/2T_Disk/Dropbox (BBSR)/Zhen_Gao/Saurav-DMR-comparison/files-for-combP")
dir.create("Comb-P-Folders", showWarnings = F)

file_names_ended_with_csv <- list.files(pattern = "\\.csv$") #dir() list all files
#file_names <- sapply(strsplit(file_names_ended_with_csv,".csv"), `[`, 1) 

for (name in file_names_ended_with_csv) {
    #name = "linear_model_result_for_mu0.025_repetition_1.csv"
    #print(name)   
    file= strsplit(name,".csv")[[1]]
    dir.create(paste0("Comb-P-Folders", "/", file), showWarnings = F) 
    current_dir = paste0("Comb-P-Folders", "/", file) 
    csv_file <- read.csv(name, header=T)
    
    order_csv_file <- csv_file[order(csv_file[,1],csv_file[,2], decreasing=F),]
    out_bed_file <- paste0("Comb-P-Folders", "/", file, "/", file, ".bed")

    write.table(order_csv_file, out_bed_file, sep="\t", quote=FALSE, row.names=FALSE)
    #system("comb-p pipelincomb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 linear_model_result_for_mu0.025_repetition_1.bed")  
    
}
