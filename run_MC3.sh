
# Project name
project_name=P4

# Get the next experiment number from the readme file
experiment_no=$((10#$(echo $(tail -n 1 readme.txt)| cut -d'-' -f 2)+1))

# Set up folders and generate the parameter files for each experiment
module load math/MATLAB/2020a
st="$(matlab -nojvm -nodisplay -nosplash -r "format compact; p = MC3_parameters(); p.set_up($experiment_no); disp(p.experiments); disp(' '); disp(p.number_of_nodes); quit")"

E="$(echo $st | rev | cut -d' ' -f3 | rev)" # Number of experiments to run
A="$(echo $st | rev | cut -d' ' -f2 | rev)" # Number of parallel jobs
F=$[experiment_no+E-1] # Number of last experiment in array

# Create entry in log 
echo "$project_name $(printf "%05d-%05d" $(($experiment_no)) $(($F)))" >> readme.txt;


#for n in $(eval echo "{$experiment_no..$F}"); do echo "$(printf "%05d\n" $(($n))) $project_name" >> readme.txt; done

# Submit batch jobs and concatenate results after job is finished
for i in $(eval echo "{$experiment_no..$F}"); do
job_id=$(sbatch --job-name=$project_name --output=Logs/$project_name.log --error=Errors/$project_name.log "--array=1-$A" job.sh $experiment_no $F $i)
echo Job $i $project_name submitted;
job_id="$(cut -d' ' -f 4 <<< $job_id)"
sbatch --dependency=afterok:$job_id concatenate.sh $experiment_no $F $i
done
