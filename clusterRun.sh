#!/bin/bash

# to submit job in parallel with one processor:
# >> LLsub clusterRun.sh 1
# to submit job in parallel with four processors:
# >> LLsub clusterRUn.sh 4


#script = "$1"
#Np = "$2"

#echo "Number of processors: $Np"
#echo "For script: $script"

matlab -nodisplay -r "eval(pRUN('testWaves_cluster',1,'grid'))"
