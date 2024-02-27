#!/bin/bash

# to submit job in parallel with one processor:
# >> LLsub clusterRun.sh 1
# to submit job in parallel with four processors:
# >> LLsub clusterRUn.sh 4


Np="$1"

echo "Number of processors: $Np"


matlab -nodisplay -r "eval(pRUN('testWaves_cluster',$Np,'grid'))"

