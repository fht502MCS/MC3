#!/bin/bash
#SBATCH --job-name=MCS_analysis
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=10:00:00
#SBATCH --error=Logs/MCS_%A.log
#SBATCH --account=chem-mcsltt-2019
#SBATCH --array=1-1

project_name=$1
module load math/MATLAB/2020a
matlab -nojvm -nodisplay -nosplash -r "abs_by_layer('$project_name');exit"


