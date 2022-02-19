#!/bin/bash
#SBATCH --job-name=MCS
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=10:00:00
#SBATCH --error=Logs/MCS_%A.log
#SBATCH --account=chem-mcsltt-2019
#SBATCH --array=1-1


experiment_no=$3 
echo Running job $3 in range $1 to $2
echo Array no $SLURM_ARRAY_TASK_ID
module load math/MATLAB/2020a
matlab -nojvm -nodisplay -nosplash -r "MC3_wrapper($experiment_no,$SLURM_ARRAY_TASK_ID);exit"


