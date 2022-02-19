
% ----- Extract results from grouped runs -----
% Recieves the N names of projects as a string
% and returns their experiment numbers as an array

function [E] = MCS_E(N)

% Find experiment numbers in each project
f = fopen('readme.txt'); reading = 1; l = 0;
C = zeros(1,length(N)); Cn = C;
while reading == 1
    line = fgetl(f);
    if ~ischar(line)
        reading = 0;
    else
        L = strsplit(line,' ');
        if length(L) == 2 
	if strcmp(L(1),N)
            Range = strsplit(string(L{2}),'-');
            C = [str2double(Range(1)) str2double(Range(2))];
        end
        end
    end
end
E = C(1):1:C(2);
end
