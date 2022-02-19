classdef photon
    properties
        number = 1; % Photon's assigned identity
        step_number = 0; % Total number of steps by photon
        step_size = 0; % Dimensionless distance until the next event
        coordinates = [0 0 0]; % Cartesian x,y,z coordinates
        direction = [0 0 1]; % Directional cosines of photon in x,y,z
        event_number = 0; % The number of absorption and scattering events
        weight = 1; % Initially assigned 1, reduced during absorption
        layer_number = 0; % Current layer number
        path = []; % Cumulative path travelled in mm
    end
    properties (Dependent)
        active; % 0/1/2 if photon is stopped/propagating/escaped
    end
    methods
        function activity = get.active(properties)
            global max_steps
            global max_events
            global max_radius
            global max_depth
            global max_layer_number
            if properties.step_number > max_steps-1
                activity = 0; % Photon has reached the maximum step number
            elseif properties.event_number > max_events-1
                activity = 0; % Photon has reached the maximum event number
            elseif sqrt((properties.coordinates(1))^2 + (properties.coordinates(2))^2) > max_radius
                activity = 0; % Photon is too far away in the X,Y-plane 
            elseif properties.weight <= 0 
                activity = 0; % At zero weight the photon is stopped
            elseif (properties.coordinates(3) <= 0 && sign(properties.direction(3))==-1) || (properties.coordinates(3) >= max_depth) 
                activity = 2; % Photon has escaped the tissue
            elseif properties.layer_number == 0 || properties.layer_number > max_layer_number 
                if properties.step_number > 0
                    activity = 2; % Photon has escaped the tissue
                else
                    activity = 1;
                end
            else
                activity = 1;
            end
        end
    end
end
