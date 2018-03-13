#!/bin/bash

for (( i = 1; i < 51; i++ )); 
#for (( i = 1; i < 2; i++ )); 
do
  cat >  ~/R/ZG-code/meth/Run_Rscript.sh <<EOF
  #LSBATCH: User input       
  #!/bin/bash
  #BSUB -P bbc
  #BSUB -J Aclust  
  #BSUB -o %J.Ac.log  
  #BSUB -e %J.Ac.err 
  #BSUB -W 24:00 
  #BSUB -n 1  
  #BSUB -q general
  #BSUB -R 'rusage[mem=4096] span[hosts=1]'             
  module load R/3.3.1
  Rscript /nethome/zxg161/For-Lily/Diff_methylation_site_02-08-2.R $i
  
EOF

  bsub < ~/R/ZG-code/meth/Run_Rscript.sh
done

