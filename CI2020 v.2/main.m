addpath('Codes/');

kernel_function = {'M15'}; % SE, M25, M15, OU
dist_type = {'Normal','T'}; % Normal, T

for m = 1:length(kernel_function)
    for n = 1:length(dist_type)
        results = GPST(kernel_function{m},dist_type{n});
    end
end