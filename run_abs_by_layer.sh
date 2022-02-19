
# Project name
project_name=P4

module load math/MATLAB/2020a

# Submit batch jobs and concatenate results after job is finished
sbatch --job-name=abs_by_layer.$project_name --output=Logs/abs_by_layer.$project_name.log --error=Errors/abs_by_layer.$project_name.log "--array=1" job_abs_by_layer.sh $project_name
echo Job abs_by_layer $project_name submitted;
