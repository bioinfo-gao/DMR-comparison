#!/bin/bash

cd "Comb-P-Folders"

FOLDERS=*

for f in $FOLDERS
do
    echo "Processing $f folder"
    cd $f
    #comb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 *.bed & PID=$!
    comb-p pipeline -c 4 --seed 1e-3 --dist 200 -p Train --region-filter-p 0.1 --annotate hg19 *.bed & PID=$!  # run in background
    sleep 2m  # sleep for 2 minute
    kill -HUP $PID  # kill the background thread

    cd ..
done


#training_input_file='No_quote_Order_Training.bed'
#comb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 No_quote_Ordered_Training.bed 
    
   
