#!/bin/bash


ave_train_runtime=0
ave_test_runtime=0

for i in {1..50}
do
    number=`echo $i`
    Train_dir='Train_'$number
    Test_dir='Test_'$number

    mkdir $Train_dir
    mkdir $Test_dir
 
    training_input_file='No_quoto_Training_'$number'_.bed'
    testing_input_file='No_quoto_Testing_'$number'_.bed'
    cd $Train_dir
        train_start=`date +%s%N` #get time in nanosecond
        comb-p pipeline -c 4 --seed 0.1 --dist 200 -p Train --region-filter-p 0.1 --annotate hg19 /media/2T_Disk/DMS/comb-p-and-data1/$training_input_file  
        train_end=`date +%s%N`
        train_runtime=$((train_end-train_start))
        train_runtime_milli=$((train_runtime/1000000))
    cd ..
    echo $train_runtime_milli >> allTrain_out.txt 2>&1 # in millisecond
    ave_train_runtime=$((ave_train_runtime+train_runtime_milli))
    
   
    cd $Test_dir
        test_start=`date +%s%N` #get time in nanosecond
        comb-p pipeline -c 4 --seed 0.1 --dist 200 -p Test --region-filter-p 0.1 --annotate hg19 /media/2T_Disk/DMS/comb-p-and-data1/$testing_input_file 
        test_end=`date +%s%N`
        test_runtime=$((test_end-test_start))
        test_runtime_milli=$((test_runtime/1000000))
    cd ..    
    echo $test_runtime_milli >> alltest_out.txt 2>&1 # in millisecond
    ave_test_runtime=$((ave_test_runtime+test_runtime_milli))

done

echo $ave_train_runtime >> ave_Train_out.txt 2>&1 # in millisecond
echo $ave_test_runtime >> ave_test_out.txt 2>&1 # in millisecond