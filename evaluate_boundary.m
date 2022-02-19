function [next_layer, distance, normal, intersect] = evaluate_boundary(p,T,k)

        % Ray intersects between photon path and layer boundaries
        D = cell(size(T,2),1); % Distances
        V = D; UV = D; 
        V1 = D; V2 = D; V3 = D; % Vectors defining triangle vertices
        for j = 1:size(T,2)
        [t, uv, v1, v2, v3] = find_intersect(p.coordinates,p.direction',...
            [T{j}(k{j}(:,1),:)],[T{j}(k{j}(:,2),:)],[T{j}(k{j}(:,3),:)]);
        V = [V;linspace(j,j,length(t))']; %UV = [UV;uv]; 
        D{j} = t; % Distances, as {layer_index}(triangle_index)
        UV{j} = uv;
        V1{j} = [V1{j};v1]; V2{j} = [V2{j};v2]; V3{j} = [V3{j};v3]; % Vertices
        end
        
        % Find closest triangle
        distance = 10000; triangle_index = 0; layer_index = 0;
        for j = 1:size(T,2)
            for i = 1:length(D{j})
                if D{j}(i) < distance && D{j}(i) > 0.0000001
                    distance = D{j}(i); 
                    triangle_index = i; 
                    layer_index = j;
                end
            end
        end
        if distance == 10000
            % Left tissue
            intersect = [0,0,0];
            normal = [0,0,0];
        else
        
        v1 = V1{layer_index}(triangle_index,:);
        v2 = V2{layer_index}(triangle_index,:);
        v3 = V3{layer_index}(triangle_index,:);
            
        % Find normal on v1,v2,v3
        try cross_p = cross(v2-v1,v3-v1); %0 = cross_p*[n-v1']
        catch
            pause(1); disp('Error in cross product');
        end
        normal = cross_p./(sqrt(cross_p(1)^2+cross_p(2)^2+cross_p(3)^2));
        
        % Barycentric coordinates to cartesian; w = 1-u-v so
            uv = UV{layer_index}(triangle_index,:);
            uvw = [1-uv(1)-uv(2),uv(1),uv(2)];
            intersect = uvw(1)*v1 + uvw(2)*v2 + uvw(3)*v3;
            
        % Determine which face of the triangle the photon is approaching
        if acosd(dot(p.direction,normal)) < 90 % Leaving current volume
            normal = normal*(-1); % Other face of triangle
            
            % Find next closest triangle
            next_distance = 10000; next_layer_index = 0;
            for j = 1:size(T,2)
                for i = 1:length(D{j})
                    if D{j}(i) < next_distance && D{j}(i) > distance && ~(j==layer_index)
                        next_distance = D{j}(i); 
                        %next_triangle_index = i; 
                        next_layer_index = j;
                    end
                end
            end
            layer_index = next_layer_index; 
        end
        end
        next_layer = layer_index;
end