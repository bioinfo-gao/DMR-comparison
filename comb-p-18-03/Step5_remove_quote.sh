#!/bin/bash
#!/bin/bash

training_pvals_file='Ordered_Training.bed' 
training_pvals_Format_file='No_quote_Ordered_Training.bed'
#awk -F',' -v OFS=',' '{gsub(/"/, "", $2); print "chr"$2,$3,$4,$5}' $training_pvals_file  > $training_pvals_Format_file
#sed -i '1s/chrchrom/chrom/' $training_pvals_Format_file

#training_pvals_tsv_file='pvals_PFC_'$number'_training-TSV.bed'
#training_pvals_sorted_file='pvals_PFC_'$number'_training-sorted.bed'
  
cat $training_pvals_file | sed 's/\"//g;s/,/\t/g' > $training_pvals_Format_file
