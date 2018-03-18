setwd("2T_Disk/Dropbox/Zhen_Gao/DMR-Saurav/code") #-------- need access to this file------------- #

data.dir1 <- "/media/2T_Disk/Dropbox/Zhen_Gao/DMR-Saurav/Saurav-DMR-comparison-files-for-combP/Comb-P-Folders"  # ----- ?????----- #
result.dir <- "2T_Disk/Dropbox/Zhen_Gao/DMR-Saurav/code/comb-p-result"

miu_values <- c(0, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40) # simulation scenarios

dir_list = list.dirs(path = "/media/2T_Disk/Dropbox/Zhen_Gao/DMR-Saurav/Saurav-DMR-comparison-files-for-combP/Comb-P-Folders")

mu0_dir     = grep(pattern = "mu0_repetition", x=dir_list)
mu0.025_dir = grep(pattern = "mu0.025_repetition", x=dir_list)
mu0.05_dir  = grep(pattern = "mu0.05_repetition", x=dir_list)
mu0.1_dir   = grep(pattern = "mu0.1_repetition", x=dir_list)
mu0.15_dir  = grep(pattern = "mu0.15_repetition", x=dir_list)
mu0.2_dir   = grep(pattern = "mu0.2_repetition", x=dir_list)
mu0.3_dir   = grep(pattern = "mu0.3_repetition", x=dir_list)
mu0.4_dir   = grep(pattern = "mu0.4_repetition", x=dir_list)

mu0_dir_names      = dir_list[mu0_dir]
mu0.025_dir_names  = dir_list[mu0.025_dir]
mu0.05_dir_names   = dir_list[mu0.05_dir]
mu0.1_dir_names    = dir_list[mu0.1_dir]
mu0.15_dir_names   = dir_list[mu0.15_dir]
mu0.2_dir_names    = dir_list[mu0.2_dir]
mu0.3_dir_names    = dir_list[mu0.3_dir]
mu0.4_dir_names    = dir_list[mu0.4_dir]

get_time_sd <- function(dir_names_vector){
    
    for(i in 1:length(dir_names_vector)){
        time[i] =as.numeric(read.table(paste0(dir_names_vector[i],"/", "run_time.txt"), sep = ","))/1000
        print(time[i])   
    }
    mean(time)
    sd(time) # = sqrt(var(time)
    result=cbind(mean(time),sd(time))
}


time_0     = get_time_sd(mu0_dir_names)
time_0.025 = get_time_sd(mu0.025_dir_names)
time_0.05  = get_time_sd(mu0.05_dir_names)
time_0.1   = get_time_sd(mu0.1_dir_names)
time_0.15  = get_time_sd(mu0.15_dir_names)
time_0.2   = get_time_sd(mu0.2_dir_names)
time_0.3   = get_time_sd(mu0.3_dir_names)
time_0.4   = get_time_sd(mu0.4_dir_names)




