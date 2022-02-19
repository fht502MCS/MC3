#!/bin/bash
#SBATCH --job-name=Concatenate_results
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=00:30:00
#SBATCH --output=Logs/concatenate.log
#SBATCH --account=chem-mcsltt-2019
#SBATCH --array=1


experiment_no=$3 
echo Running job $3 in range $1 to $2
echo Array no $SLURM_ARRAY_TASK_ID
module load math/MATLAB/2020a
matlab -nojvm -nodisplay -nosplash -r "MCS_concatenate($experiment_no);exit"


