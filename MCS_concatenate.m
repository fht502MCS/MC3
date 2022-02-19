% ----- Script to concatenate MCS data from multiple runs -----

function MCS_concatenate(E)

folder_name = strcat('Data/Exp',sprintf('%05d',E),'/'); save_to = strcat('Exp',num2str(E),'.mat');

parts = strsplit(folder_name,'/');
save_to = parts{end-1};


%file_names = dir('Data/Exp1/*.mat');
dir_name = strcat(folder_name,'*.mat');
file_names = dir(dir_name); disp(file_names); 
file_char = char({file_names.name}); disp(file_char);
file_num = ismember(file_char(:,1),'0123456789');
file_names = file_char(file_num,:);
file_no = size(file_names,1);

for i = 1:file_no
    a = file_names(i,:);
    disp(strcat('Processing file:',{' '},string(a)));
    a = strcat(folder_name,a);
    b = load(a);
    if i == 1
        simulation = b.simulation;
        tissue = b.tissue;
        boundaries = b.boundaries;
        rdr = b.rdr;
	psf = b.psf;
	Td = b.Td;
	coordinates = b.coordinates;
        paths = b.paths;
	all_paths = b.all_paths;
	R_unscat = b.R_unscat;
	T_unscat = b.T_unscat;
	%lost = b.lost;
	edges = b.edges;
	idx = b.idx;
	bin_size = b.bin_size;
	R_layers = b.R_layers;
	T_layers = b.T_layers;
	Escaped_bounds = b.Escaped_bounds;
	abs_coords = b.abs_coords;
	abs_weight = b.abs_weight;
    abs_layer = b.abs_layer;
    abs_id = b.abs_id;
    else
        simulation.number_of_photons = simulation.number_of_photons+b.simulation.number_of_photons;
        %simulation.wavelength = [simulation.wavelength,[b.simulation.wavelength]];
        %tissue = [tissue,[b.tissue]];
        %boundaries = [boundaries,[b.boundaries]];
	rdr = [rdr;[b.rdr(2,:)]];
	psf = cat(3,psf,b.psf(:,:,2));
        coordinates = [coordinates,[b.coordinates]];
        paths = [paths;b.paths];
	all_paths = [all_paths;b.all_paths];
	Td = cat(3,Td,b.Td(:,:,2));
	R_unscat = R_unscat + b.R_unscat;
	T_unscat = T_unscat + b.T_unscat;
	R_layers = [R_layers;[b.R_layers]];
	T_layers = [T_layers;[b.T_layers]];
	Escaped_bounds.weight = [Escaped_bounds.weight;b.Escaped_bounds.weight];
	Escaped_bounds.coordinates = [Escaped_bounds.coordinates;b.Escaped_bounds.coordinates];
	Escaped_bounds.r = [Escaped_bounds.r;b.Escaped_bounds.r];
	abs_coords = [abs_coords;b.abs_coords];
	abs_weight = [abs_weight;b.abs_weight];
    abs_layer = [abs_layer;b.abs_layer];
	abs_id = [abs_id;b.abs_id];
	%lost = lost + b.lost;
    end
end

save_file = strcat(folder_name,save_to);
save(save_file,'simulation','tissue','boundaries','rdr','paths',...
    'coordinates','idx','bin_size','edges','psf','Td','R_unscat',...
    'T_unscat','R_layers','T_layers','Escaped_bounds','abs_coords',...
    'abs_weight','abs_layer','abs_id');
disp('Files concatenated');
end
