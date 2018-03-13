#!/bin/bash
    # for training data format processing
    training_pvals_file='pvals_PFC_training.bed' 
    training_pvals_Format_file='pvals_PFC_training-Format.bed'
    awk -F',' -v OFS=',' '{gsub(/"/, "", $2); print "chr"$2,$3,$4,$5}' $training_pvals_file  > $training_pvals_Format_file
    sed -i '1s/chrchrom/chrom/' $training_pvals_Format_file

    training_pvals_tsv_file='pvals_PFC_TSV.bed'
  
    cat $training_pvals_Format_file | sed 's/\"//g;s/,/\t/g' > $training_pvals_tsv_file
    
    
