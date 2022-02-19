classdef MC3_parameters
    properties 
        number_of_photons = 1000; % Photons used in each node
        number_of_nodes = 50; % Number of parallel processes
        method = 'element'; % Select 'element' or 'array': 1 experiment per 
        % index or experiments covering all combinations
        
        % Loop over properties: {value [range] value} loops over 2nd layer 
        %load('wavelength_steps');
        sweep = 'sweep_steps'; % Load from file
        wavelength = 'wavelength_steps'; %linspace(560,560,10)'}; % 1D only {[660]};
        layers = 'P4_Mesh'; % Layer thickness
        refractive_index = {linspace(1.4,1.4,131)',linspace(1.4,1.4,131)',linspace(1.4,1.4,131)',linspace(1.4,1.4,131)',linspace(1.4,1.4,131)',linspace(1.4,1.4,131)'};
        muA = 'muA_steps';
        muSr = 'muSr_steps'; 
        g = {linspace(0.88,0.88,131)',linspace(0.94,0.94,131)',linspace(0.999,0.999,131)',linspace(0.96,0.96,131)',linspace(0.87,0.87,131)',linspace(0.91,0.91,131)'};
               
        % Not looped over:
        max_steps = 10000; % Max steps before a photon is stopped
        max_events = 10000; % Max scattering or absorption events
        max_radius = 50; % Max distance from the Z-axis, mm
        max_angle = 43; % Angle in degrees in which the camera detects   
        max_depth = 75;
	threshold_weight = 0.0001; m = 10; % Weight bounds
    end
    properties (Dependent)
        experiments;
        MuA;
        MuSr;
        Wavelength;
        Layers;
        UseTriMesh;
    end
    methods
        % Find number of experiments needed to loop over properties
        function MuA = get.MuA(properties)
            if ischar(properties.muA)
                muA_list = load(properties.sweep,properties.muA);
                F = fieldnames(muA_list); MuA = {};
                for i = 1:size(muA_list.(F{1}),2)
                MuA = [MuA,muA_list.(F{1})(:,i)];
                end
            else
                MuA = properties.muA;
            end
        end
        function MuSr = get.MuSr(properties)
            if ischar(properties.muSr)
                muSr_list = load(properties.sweep,properties.muSr);
                F = fieldnames(muSr_list); MuSr = {};
                for i = 1:size(muSr_list.(F{1}),2)
                MuSr = [MuSr,muSr_list.(F{1})(:,i)];
                end
            else
                MuSr = properties.muSr;
            end
        end
        function Wavelength = get.Wavelength(properties)
            if ischar(properties.wavelength)
            wavelength_list = load(properties.sweep,properties.wavelength);
            F = fieldnames(wavelength_list);
            Wavelength = {wavelength_list.(F{1})};
            else
            Wavelength = properties.wavelength;
            end
        end
        function Layers = get.Layers(properties)
            if ischar(properties.layers)
            layers = load(properties.layers);
            Layers = layers.T;
            else
            Layers = properties.layers;
            end
        end
        function UseTriMesh = get.UseTriMesh(properties)
            if ischar(properties.layers)
            UseTriMesh = 1;
            else
            UseTriMesh = 0;
            end
        end

	function exp = get.experiments(properties)
            L = size(properties.Layers,2); 
            if size(properties.MuA,2) == L && ...
                    size(properties.MuSr,2) == L &&...
                    size(properties.refractive_index,2) == L &&...
                    size(properties.g,2) == L
                if strcmp(properties.method,'array')
                exp = length(properties.Wavelength);
                for i = 1:L
                    exp = exp*size(properties.layers{i},2)*...
                    size(properties.refractive_index{i},2)*...
                    size(properties.MuA{i},2)*...
                    size(properties.MuSr{i},2)*...
                    size(properties.g{i},2);
                end
                elseif strcmp(properties.method,'element')
                    A = size(properties.MuA{1},1); 
                if size(properties.Wavelength{1},1) == A && ...
                    size(properties.refractive_index{1},1) == A && ...
                    size(properties.MuA{1},1) == A && ...
                    size(properties.MuSr{1},1) == A && ...
                    size(properties.g{1},1) == A
                    exp = A;
                else
                    exp = 0; disp('Properties not assigned correctly');
                end
                end
            else
                exp = 0; disp('Properties not assigned correctly');
            end               
        end
        % Generate data folders and the parameter file in each
        function set_up(properties,E)
            simulation = struct(...
                'number_of_photons',properties.number_of_photons,...
                'number_of_nodes',properties.number_of_nodes,... 
                'useTriMesh',properties.UseTriMesh);
            boundaries = struct(...
                'max_steps',(properties.max_steps),...
                'max_events',(properties.max_events),...
                'max_radius', (properties.max_radius),...
                'max_angle', (properties.max_angle),...
                'max_depth', (properties.max_depth),...
                'threshold_weight', (properties.threshold_weight),...
                'm', (properties.m));
if strcmp(properties.method,'array')
            R = {properties.refractive_index{:},...
                properties.MuA{:},...
                properties.MuSr{:},...
                properties.g{:},...
                properties.layers{:}};
            
	    elements = R; %cell array with N vectors to combine
            combinations = cell(1, numel(elements));
            [combinations{:}] = ndgrid(elements{:});
            combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); 
            result = [combinations{:}]; 
            L = size(properties.layers,2);
end
            for w = 1:length(properties.Wavelength{:})
            for i = 1:size(properties.Layers,2)
                 if strcmp(properties.method,'array')
                     mkdir(strcat('Data/Exp',sprintf('%05d',(E-1+i+(w-1)*...
                     (properties.experiments/length(properties.Wavelength{:}))))));
                    tissue = struct(...
                        'wavelength',(properties.Wavelength(w)),...
                        'refractive_index',(result(i,1:L)),...
                        'muA',(result(i,(L+1):(2*L))),...
                        'muSr',(result(i,(2*L+1):(3*L))),...
                        'g',(result(i,(3*L+1):(4*L))),...
                        'layers',(result(i,(4*L+1):(5*L))));
                elseif strcmp(properties.method,'element')
                    if i == 1
                        mkdir(strcat('Data/Exp',sprintf('%05d',(E-1+w))));
                    tissue = struct(...
                        'wavelength',(properties.Wavelength{i}(w)),...
                        'refractive_index',(properties.refractive_index{i}(w)),...
                        'muA',(properties.MuA{i}(w)),...
                        'muSr',(properties.MuSr{i}(w)),...
                        'g',(properties.g{i}(w)),...
            			'layers',({properties.Layers}));%(properties.Layers{i}));
                    else
                        
                        tissue.refractive_index = [tissue.refractive_index,properties.refractive_index{i}(w)];
                        tissue.muA = [tissue.muA,properties.MuA{i}(w)];
                        tissue.muSr = [tissue.muSr,properties.MuSr{i}(w)];
                        tissue.g = [tissue.g,properties.g{i}(w)];
                        %tissue.layers = properties.Layers;%[properties.Layers{i}];
                        %tissue = [tissue,t];
                    end
                 end
            end
                parameters = struct();
                parameters.simulation = simulation;
                parameters.tissue = tissue;
%parameters.tissue.layers = properties.Layers;
                parameters.boundaries = boundaries;
                if strcmp(properties.method,'array')
                filename = strcat('Data/Exp',sprintf('%05d',(E-1+i+(w-1)*...
                     round(properties.experiments/length(properties.Wavelength{:})))),'/parameters');
                else
                    filename = strcat('Data/Exp',sprintf('%05d',(E-1+w)),'/parameters');
                end
                save(filename,'parameters');
            end
	    end
        
    end
end
