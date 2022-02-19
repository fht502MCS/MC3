
%        Monte Carlo simulation of light transport through tissue         

function [coordinate_store, path_store, Rdr, Rd, delta_r, idx, edges,...
    Td, R_unscat, T_unscat, Tdr, R_layers, T_layers, Escaped_bounds,...
    Roulette_weight, R_e, T_e, all_paths, abs_coords, abs_weight,...
    abs_layer, abs_id] = MC3_function(simulation, tissue, boundaries)

% Make some loaded boundaries global for passing to the photon class
global max_steps, max_steps = boundaries.max_steps;
global max_events, max_events = boundaries.max_events;
global max_radius, max_radius = boundaries.max_radius;
global max_depth, max_depth = boundaries.max_depth;
global max_layer_number, max_layer_number = size(tissue.layers,2);

% ----- Set more parameters -----
tissue.muT = tissue.muA + tissue.muSr;
clear_layer = 0; n1 = 1; % Refractive index of air %1.00027717; 
n3 = 1; % Set n3 for ambient medium

if (simulation.useTriMesh == 1)
    k = {}; T = tissue.layers;
    for i = 1:size(T,2)
        k{i}=boundary(T{i});
    end
else
    layer_edges = zeros(1,length(tissue.layers)+1);
    for n = 1:length(tissue.layers)
        layer_edges(n+1) = sum(tissue.layers(1:n));
    end
end

% ----- Storage allocation -----
path_store = []; coordinate_store = {}; R_unscat = 0; T_unscat = 0;
all_paths = []; abs_coords = {}; abs_weight = {}; abs_layer = {}; abs_id = {};

% Radial diffuse reflectance % Cuccia: radial bins spacing 0.09 mm
Nr = ceil(max_radius/0.09); delta_r = max_radius/Nr; % Capture radius/Bin number 
Rdr = zeros(2,Nr); idx = zeros(1,Nr);
for i = 0:Nr-1   % Wang
    idx(1,i+1) = ((i+0.5)+(1/(12*(i+0.5))))*delta_r; % indices
end
for i = 0:Nr
    r = ((i)+(1/(12*(i))))*delta_r; % edges
    Rdr(1,i+1) = r;
end
Rdr(1,1) = 0; % when using edges
Tdr = Rdr;
edges = Rdr(1,:); % Rdr(1,1:end-1) = layer_index;

% PSF of diffuse reflectance
N = ceil(max_radius/0.09); b = linspace(0,max_radius,N); 
b = [fliplr(-b),b(2:end)]; Rd = zeros(length(b),length(b),2);
for i = 1:length(b)
   gx = b(i);
   for j = 1:length(b)
       gy = b(j); gr = sqrt(gx^2 + gy^2); Rd(i,j,1) = gr;
   end
end
Td = Rd;

R_e = struct(); R_e.c = [0 0 0]; R_e.p = [0]; T_e = R_e;

R_layers = zeros(1,max_layer_number); % Internal reflectance in layer
T_layers = R_layers; % Transmittance from each layer
Escaped_bounds = {}; % Exceeded max radius
Escaped_bounds.coordinates = [0 0 0];
Escaped_bounds.weight = [0];
Escaped_bounds.r = [0];
Roulette_weight = 0; % Weight added to the n*1 photon weight

% ----- Run the procedure -----
rng('shuffle'); % Shuffle random numbers
%h = waitbar(0,'Running simulation...');
for n = 1:simulation.number_of_photons

    % ----- Initialise a new photon -----
    p = photon(); p.number = n;
    photon_store = [p.coordinates]; % Log initial position
    p.path = zeros(1,max_layer_number); % Record path in each layer    
    
    % ----- Specular reflectance -----
    if (simulation.useTriMesh == 1) % Launch virtual photon to determine start
        p0 = photon(); p0.coordinates = [0,0,0]; p0.step_size = 2;
        [next_layer, ~, ~, intersect] = evaluate_boundary(p0,T,k);
        n2 = tissue.refractive_index(next_layer);
    else
        intersect = [0,0,0];
        n2 = tissue.refractive_index(1); 
        next_layer = 1;
    end
    if clear_layer == 0 % TODO: Include medium in layers with alt. method
        Rsp = (n1-n2)^2/(n1+n2)^2;                
    else % If the photon travels through a clear medium   
        r1 = (n1-n2)^2/(n1+n2)^2; r2 = (n3-n2)^2/(n3+n2)^2;
        Rsp = r1 + ((1-r1)^2*r2)/(1-r1*r2);
    end            
    % Specular reflectance decreases the weight      
    p.weight = p.weight - Rsp; 
    
    % First step moves photon to boundary of first layer
    p.coordinates = intersect;
    p.layer_number = next_layer;    
    
    while(p.active == 1) 
        if p.step_size == 0 % Reset step size for new or scattered photons
            p.step_size = -log(rand);        
        end     
        
        % ----- Layer boundaries -----
        % A photon can be internally reflected or transmit across boundary
        if (simulation.useTriMesh == 1)
            [next_layer, distance, normal, intersect] = evaluate_boundary(p,T,k);
            c = dot(p.direction,-normal);
        else
            I = layer_edges(p.layer_number); J = layer_edges(p.layer_number + 1); 
            % Distance from photon location to the current layer boundaries I,J
            if p.direction(3) < 0
                distance = (I-p.coordinates(3))/p.direction(3);
            elseif p.direction(3) == 0
                distance = inf;
            elseif p.direction(3) > 0
                distance = (J-p.coordinates(3))/p.direction(3);
            end
            intersect = p.coordinates + p.direction*distance;
            next_layer = p.layer_number+sign(p.direction(3));
            c = abs(p.direction(3));
        end
        
        if p.active == 1
        % Check if the dimensionless stepsize is greater than the distance
        interaction = distance*tissue.muT(p.layer_number) <= p.step_size;
        if interaction % If true, the photon hits the tissue boundary 
            
            % Photon is moved to the boundary and the step size updated
            p.coordinates = intersect;          
            p.step_number = p.step_number + 1;
            p.step_size = p.step_size - distance*tissue.muT(p.layer_number);
            p.path(p.layer_number) = p.path(p.layer_number) + distance;                 
            if p.active ~= 0 % If still within boundaries

            % ----- Internal reflection vs transmittance -----
            alpha_i = acosd(c);
            nI = tissue.refractive_index(p.layer_number); % Current layer
            if next_layer == 0 || next_layer > max_layer_number % Leaving tissue
                nT = 1; % Refractive index of air               
            else 
                nT = tissue.refractive_index(next_layer);
            end            
            alpha_t = asind((nI*sind(alpha_i))/nT); % Snell's law 
            if alpha_i > asind(nT/nI)
                internal_reflectance = 1;
            else % Fresnel's formulas:
                internal_reflectance = (1/2)*((sind(alpha_i-alpha_t))^2/...
                ((sind(alpha_i+alpha_t))^2)+((tand(alpha_i-alpha_t))^2)/...
                ((tand(alpha_i+alpha_t))^2));
            end          
            
            % Decide if the photon is internally reflected or transmitted
            if rand <= internal_reflectance % Internally reflected
                if (simulation.useTriMesh == 1)
                    p.direction = p.direction-2*(dot(p.direction,normal))*normal; % Direction after reflection
                else
                    p.direction = p.direction.*[1 1 -1]; % Reverse z-direction
                end
                R_layers(p.layer_number) = R_layers(p.layer_number) + p.weight;
            else % Transmitted:
                T_layers(p.layer_number) = T_layers(p.layer_number) + p.weight;
                if next_layer == 0 || next_layer > max_layer_number % Photon escaped the tissue
                    nT = 1; % Refractive index of air
                else 
                    nT = tissue.refractive_index(next_layer); % To 
                end 
                % Find the new directions
                nI = tissue.refractive_index(p.layer_number); % From
                if (simulation.useTriMesh == 1)
                    r = nI/nT;
                    v_refract = r*p.direction + (r*c-sqrt(1-(r^2)*(1-c^2)))*normal;
                    p.direction = v_refract;               
                else
                    x_dir = p.direction(1)*nI/nT;
                    y_dir = p.direction(2)*nI/nT;
                    z_dir = sign(p.direction(3))*cosd(alpha_t);   
                    p.direction = [x_dir, y_dir, z_dir]; 
                end
                p.layer_number = next_layer; 
            end 
            end
        else % If not hitting the layer boundary, the photon moves by s/muT
            % ----- Photon movement -----   
            p.coordinates = p.coordinates + (p.direction)*p.step_size/...
                tissue.muT(p.layer_number);
            p.path(p.layer_number) = p.path(p.layer_number) + ...
                p.step_size/tissue.muT(p.layer_number); 
            p.step_number = p.step_number + 1; p.step_size = 0; % Reset
            if p.active ~= 0 % If still within boundaries

            % ----- Photon absorption -----      
            weight_change = (tissue.muA(p.layer_number)/...
                tissue.muT(p.layer_number))*p.weight;
            p.weight = p.weight - weight_change;
            abs_coords = [abs_coords;p.coordinates];
            abs_weight = [abs_weight;weight_change];
            abs_layer = [abs_layer;p.layer_number];
            abs_id = [abs_id;p.number];

            % ----- Photon scattering -----                                    
            if tissue.g == 0 % Isotropic scattering case
                cos_theta = 2*rand -1;
            else % Anisotropic case, Henyey-Greenstein phase function
                g = tissue.g(p.layer_number);
                cos_theta = (1/(2*g))*(1+g^2-((1-g^2)/(1-g+2*g*rand))^2);
            end
            theta = acos(cos_theta); % The deflection angle
            phi = 2*pi*rand; % The azimuthal angle is uniformly distributed 
            % Calculate the new direction based on these angles: [1]
            if abs(p.direction(3)) < 0.99999 
                x_dir = sin(theta)*(p.direction(1)*p.direction(3)*cos(phi)-...
                    p.direction(2)*sin(phi))/(sqrt(1-p.direction(3)^2))+...
                    p.direction(1)*cos(theta);
                y_dir = sin(theta)*(p.direction(2)*p.direction(3)*cos(phi)+...
                    p.direction(1)*sin(phi))/(sqrt(1-p.direction(3)^2))+...
                    p.direction(2)*cos(theta);
                z_dir = -sin(theta)*cos(phi)*sqrt(1-p.direction(3)^2)+...
                    p.direction(3)*cos(theta);
            else % Photon direction is close to the z-axis [Wang 3.25] 
                x_dir = sin(theta)*cos(phi);                                
                y_dir = sin(theta)*sin(phi);
                z_dir = sign(p.direction(3))*cos(theta);
            end
            p.direction =  [x_dir y_dir z_dir];
            end
        end
        end
        if p.active == 1 % Russian roulette
            if p.weight < boundaries.threshold_weight % Small weight 
                if rand <= 1/boundaries.m
                    Roulette_weight = Roulette_weight + (boundaries.m*p.weight - p.weight);
                    p.weight = boundaries.m*p.weight; 
                else
                    Roulette_weight = Roulette_weight + p.weight;
                    p.weight = 0;
                end
            end
        end            
        % Store the coordinates for each photon step
        photon_store = [photon_store ; p.coordinates];   
        
    end % Photon no longer propagating
    coordinate_store{end+1} = photon_store;
    
    % --- Update with results from the last photon ---
    all_paths = [all_paths ; p.path]; % From all photons
    if p.active == 2 % Photon escaped the tissue        
        % Store the path lengths of remitted photons
        if p.layer_number == 0
            path_store = [path_store ; p.path]; 
        end	
        % Score escaped photon weight into coordinate array
        Escaped_bounds.coordinates = [Escaped_bounds.coordinates;p.coordinates];
        Escaped_bounds.weight = [Escaped_bounds.weight;p.weight];
        Escaped_bounds.r = [Escaped_bounds.r;(-1)];      
    end
 %   disp([num2str(n),'/',num2str(simulation.number_of_photons)]);
%    waitbar(n/simulation.number_of_photons,h)
    
end
%close(h); 
end
