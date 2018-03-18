# Input: 
# Combp.ranges = a dataframe of Combp results
# aclustfile = results from A-cluater

processCombpResults <- function(combp.ranges, aclustfile, mval_each, r, cpg_number_cutoffval) ##newly added last two arguments
{
    if(is.null(combp.ranges)==FALSE){
        
        # save result with significant DMRs that contains length greater than four is taken
        results.ranges11_save_for_Corr_analysis <- combp.ranges[which(combp.ranges$z_p < 0.05 & combp.ranges$n_probes > cpg_number_cutoffval),]  ##newly added for Correlation analysis
        write.table(results.ranges11_save_for_Corr_analysis,file=paste("combp_result_final_for_miu",mval_each,"_rep",r,".csv",sep=""),sep=",",row.names=F,col.names=TRUE) ##newly added for Correlation analysis
        
        # make ranges for Combp results (newly formatted)
        require(IRanges)
        require(GenomicRanges)
        combp.ranges.Iranges <- IRanges(combp.ranges$start, combp.ranges$end)##as class(combp.ranges$start) and class(combp.ranges$end) are both in "integer" class
        query <- GRanges(seqname = combp.ranges$`#chrom`, ranges = combp.ranges.Iranges)#Combp results
        
        ###take only region or cluster start and end info eliminating cpgs name, coordinate and beta values for using overlap function
        ##reformatted the first line
        res_aclcpg_repeach_region <- cbind(aclustfile$Clusternumber,
                                           as.character(levels(aclustfile$chromosome)[aclustfile$chromosome]),
                                           aclustfile$start_position,
                                           aclustfile$end_position,
                                           as.character(levels(aclustfile$actual)[aclustfile$actual]))
        
        colnames(res_aclcpg_repeach_region) <- c("Clusternumber","chromosome","start_position","end_position","actual")
        res_aclcpg_repeach_region_uniq <- as.data.frame(res_aclcpg_repeach_region[!duplicated(res_aclcpg_repeach_region[,1]),]) ##unique cluster information of 500 random clusters from Aclust results
        
        start <- as.numeric(levels(res_aclcpg_repeach_region_uniq$start_position)[res_aclcpg_repeach_region_uniq$start_position])
        end <- as.numeric(levels(res_aclcpg_repeach_region_uniq$end_position)[res_aclcpg_repeach_region_uniq$end_position])
        aclustfile.Iranges <- IRanges(start, end)
        subject <- GRanges(seqname=res_aclcpg_repeach_region_uniq$chromosome, ranges=aclustfile.Iranges)#aclustfile (simAclust.file) results
        
        # overlap
        overlapinfo <-  as.data.frame(findOverlaps(query, subject, maxgap=0L, minoverlap=1L, type="any", select="all") )
        overlapinfo$dmr.order <- as.numeric(overlapinfo$queryHits)
        overlapinfo$aclust.order <- as.numeric(overlapinfo$subjectHits)
        
        # merge with combp info
        results <- combp.ranges
        
        results$dmr.row <- as.numeric(seq.int(nrow(results)))
        results$predicted <- ifelse(results$z_p < 0.05, "positive", "negative")
        
        dmrs.overlap <- merge(x=overlapinfo, y=results, by.x="dmr.order", by.y="dmr.row", all=TRUE, all.x=TRUE, all.y=TRUE)
        
        # merge with aclust info
        aclust <- res_aclcpg_repeach_region_uniq
        
        aclust$aclust.row <- as.numeric(seq.int(nrow(aclust)))
        
        dmrs.aclust.overlap <- merge(x=dmrs.overlap, y=aclust, by.x="aclust.order", by.y="aclust.row", all=TRUE)
        
        # classify into different status
        all <- dmrs.aclust.overlap
    }
    
    if(is.null(combp.ranges)==TRUE){
        all <- aclustfile
        all$predicted <- "negative"
    }
    
    
    all$predicted [is.na(all$predicted)] ="negative"
    all$actual [is.na(all$actual)] = "negative"
    
    all$status[all$actual=="positive" & all$predicted=="positive"] = "TP"
    all$status[all$actual=="positive" & all$predicted=="negative"] = "FN"
    
    all$status[all$actual=="negative" & all$predicted=="positive"] = "FP"
    all$status[all$actual=="negative" & all$predicted=="negative"] = "TN"
    
    # save result ##newly added
    
    write.table(all,file=paste("mergedresult_of_Combp_and_aclust_for_miu",mval_each,"_repition",r,".csv",sep=""),sep=",",row.names=F,col.names=T)##newly added
    
    # Frequency count of each status according to Combp  ##newly added
    all_aclstno_status<-all[,c("Clusternumber","status")]
    nrow(all_aclstno_status) ##3065
    all_aclstno_status_uniq<-all_aclstno_status[!duplicated(all_aclstno_status[c(1,2)]),]## duplicacy removal of aclust clusternumber and status combinedly
    nrow(all_aclstno_status_uniq) ##3063
    
    status.count <- as.data.frame(table(all_aclstno_status_uniq$status))
    
    # power and pecision compute ##newly added
    require(tidyr)
    temp.one <- spread(data=status.count, key=Var1, value=Freq)
    if(length(grep("FP", colnames(temp.one)))==0) {temp.one$FP = 0}
    if(length(grep("TP", colnames(temp.one)))==0) {temp.one$TP = 0}
    if(length(grep("TN", colnames(temp.one)))==0) {temp.one$TN = 0}
    if(length(grep("FN", colnames(temp.one)))==0) {temp.one$FN = 0}
    
    temp.one$Power <- temp.one$TP/(temp.one$TP + temp.one$FN)
    temp.one$Precision <- temp.one$TP/(temp.one$TP + temp.one$FP)
    
    temp.one <- temp.one[, c("FN", "FP", "TN", "TP", "Power", "Precision")] #added 11-22, need to reorder
    
    return (temp.one)
}