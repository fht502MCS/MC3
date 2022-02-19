
% --- Extract absorption data by depth ---

function [] = abs_by_layer(experiment_name)

% Find experiment numbers in each project
E = MCS_E(experiment_name); 

project_name = ['E',num2str(E)];
parameter_file = strcat('Data/Exp',sprintf('%05d',E(1)),'/parameters.mat');
if ~~exist(parameter_file)
load(parameter_file); 
z_bins = linspace(0,parameters.boundaries.max_depth,101);
Z = z_bins(1:end-1);

%A = zeros(length(E),100); 
L = zeros(length(E),length(parameters.tissue.muA));
nodes = parameters.simulation.number_of_nodes; P = zeros(1,parameters.simulation.number_of_nodes);
for e = 1:length(E)
    filename = strcat('Data/Exp',sprintf('%05d',E(e)),'/',num2str(e)); 
    z = []; r = []; a = []; l = []; photons = 0; 
    for n = 1:parameters.simulation.number_of_nodes
        %disp(strcat('___',num2str(n),'/',num2str(parameters.simulation.number_of_nodes)));
        f = strcat('Data/Exp',sprintf('%05d',E(e)),'/',num2str(n),'.mat');
        try load(f,'abs_coords');
            C = cell2mat(abs_coords);
	    z = C(:,3); z(z>50) = 0; z(z>0) = 1;
	    r = (sqrt(C(:,1).^2+C(:,2).^2)); % Radius from photon launch axis
	    r(r>50) = 0; r(r>0) = 1; 
            %R = reshape(R,3,length(R)/3);
            %z = [z;R(3,:)']; % z-coordinate only
            %r = [r;(sqrt(R(1,:).^2+R(2,:).^2))']; % Radius from photon launch axis
            photons = photons + parameters.simulation.number_of_photons;
        catch
            % Skip
            nodes(1,e) = 1;
        end
        try load(f,'abs_weight');
            aa = [abs_weight{:}]; 
	    aa(z==0) = 0; aa(r==0)=0;
            a = [a;aa'];  
        catch
            % Skip
            nodes(2,e) = 1;
        end
        try load(f,'abs_layer');
            ll = [abs_layer{:}]; 
            l = [l;ll'];
        catch
        end
    end
% idx = zeros(length(z_bins)-1,length(z)); Abs = zeros(1,length(z_bins)-1);
% if ~isempty(a)
% for i = 1:length(z_bins)-1
%     idx(i,:) = (z>=z_bins(i))&(z<z_bins(i+1)); 
%     s = sum(a(logical((idx(i,:))))',2);
%     if ~isempty(s)
%     Abs(i) = s;
%     end
% end
% end

layer_abs = zeros(1,length(parameters.tissue.muA));
for i = 1:length(parameters.tissue.muA)
    layer_abs(i) = sum(a(l==i));
end
L(e,:) = layer_abs;

%A(e,:) = Abs;
P(e) = photons;
disp(strcat(num2str(e),'/',num2str(length(E))));

end
%save(strcat('Results/',project_name,'_A_Z_L'),'A','Z','L','P','nodes');
save(strcat('Results/',experiment_name,'_L_2'),'L','P','nodes');
end
end
