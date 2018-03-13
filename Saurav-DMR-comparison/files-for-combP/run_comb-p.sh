#!/bin/bash

cd "Comb-P-Folders"

FOLDERS=*

for f in $FOLDERS
do
echo "Processing $f folder"
#cd $f
cd $f
#touch 
#ls *.bed > bed_name
#grep ".bed" > sed_file
#head sed_file
#ls
comb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 *.bed
cd ..
done


#training_input_file='No_quote_Order_Training.bed'
#comb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 No_quote_Ordered_Training.bed 
    
   
