% ---------- Monte Carlo Simulations ----------

% How to run
 - Adjust parameters in MCS_parameters.m Tissue 
properties are given as lists, with each entry 
corresponding to a tissue layer.
Assign n x m property arrays, where n = no. of 
layers, m = no. of properties. Select 'element' 
or 'array': giving either m experiments (each 
with n layers) OR m_1*m_2*...m_i experiments 
covering all combinations
 - Edit the job.sh script to specify the number 
of runs required for these parameters.
The -t entry gives the number of files and the 
labels the files will get.
 - If running parallell simulations, concatenate
the results when the simulation is finished using 
MCS_concatenate. Edit the folder name and file 
name the results should be saved to before 
submitting.

% Scripts
MCS_function.m     contains the main simulation.
MCS_parameters.m   sets the parameters for the 
		   current run, calls the main 
		   fucntion and saves the data.
MCS_parameters.job submits the parameter function 
	 	   to run. Allows for parallell 
		   processing of n photons on m 
		   nodes.
photon.m 	   a class called by the main 
		   function, it contains the 
		   properties of the current 
		   photon.
MCS_concatenate.m  concatenates the outputs of
		   a simulation, m runs of n 
		   photons each provides one 
		   dataset from n*m photons.


Experiments
Initialise 00000-00000