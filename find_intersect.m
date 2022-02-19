function [t, uv, v1, v2, v3] = find_intersect(O,D,V0,V1,V2)

% O = repmat(O',size(V0,1),1);
% D = repmat(D,size(V0,1),1);
O = O;
% V0, V1 and V2 are n-x-3 coordinates
E1 = V1-V0;
E2 = V2-V0;

Ep = 0.000001; % Cut-off for determinant division


t = []; uv = []; v1 = []; v2 = []; v3 = []; 
for i = 1:size(V0,1)
    T = O-V0(i,:); % Distance from V0 to ray origin
    P = cross(D,E2(i,:));
    Q = cross(T,E1(i,:));
    det = (dot(cross(D,E2(i,:)),E1(i,:)));
    if ~(det > -Ep && det < Ep)
        inv_det = 1/det;
        U = dot(T,P)*inv_det;
        if ~(U < 0 || U > 1)
            V = dot(D,Q)*inv_det;
            if ~(V < 0 || U+V > 1)
                t = [t; dot(E2(i,:),Q)*inv_det];
                uv = [uv; [U V]];
                v1 = [v1; V0(i,:)];
                v2 = [v2; V1(i,:)];
                v3 = [v3; V2(i,:)];
            end
        end
    end
end

end
