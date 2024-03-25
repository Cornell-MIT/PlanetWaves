#!/bin/bash
#SBATCH --output=runWaves_%j.log
#SBATCH -n 3
#SBATCH -N 1

# load module (no paralell)
matlab -nodisplay -singleCompThread -r "run('$1');exit;"


