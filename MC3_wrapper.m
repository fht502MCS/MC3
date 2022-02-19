function [] = MC3_wrapper(exp_no,n)


load(strcat('Data/Exp',sprintf('%05d',exp_no),'/parameters.mat'));

simulation = parameters.simulation;
tissue = parameters.tissue;
boundaries = parameters.boundaries;

[coordinates, paths, rdr, psf, bin_size, idx, edges, Td, R_unscat, T_unscat, Tdr, R_layers, T_layers, Escaped_bounds, Roulette_weight, R_e, T_e, all_paths, abs_coords, abs_weight, abs_layer, abs_id] = MC3_function(simulation, tissue, boundaries);

% ----- Save the data and parameters -----
file_string = strcat('Data/Exp',sprintf('%05d',exp_no),'/',num2str(n));
save(file_string, 'coordinates', 'paths', 'rdr', 'psf', 'bin_size',...
    'idx', 'edges', 'Td','R_unscat','T_unscat','Tdr','R_layers','T_layers','Escaped_bounds','Roulette_weight','R_e','T_e','all_paths','abs_coords','abs_weight','abs_layer','abs_id','simulation', 'tissue', 'boundaries'); 

end
