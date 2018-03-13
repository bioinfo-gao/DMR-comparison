#!/bin/bash

#training_input_file='No_quote_Order_Training.bed'
start=`date +%s%N` #get time in nanosecond
comb-p pipeline -c 4 --seed 1e-1 --dist 500 -p Train --region-filter-p 0.1 --annotate hg19 No_quote_Ordered_Training.bed 
    
end=`date +%s%N`
runtime=$((end-start))
echo $((runtime/1000000)) # in millisecond

# 37 second