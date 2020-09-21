function [kernel_function] = saveResults(kernel_function,dist_type,results)

% Storing the results:
path = ['Outputs/',kernel_function,'_',dist_type];
if exist(path,'dir') == 7
    n = 0;
    DET = 1;
    while DET == 1
        n = n + 1;
        path = ['Outputs/',kernel_function,'_',dist_type,'(',num2str(n),')'];
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end
mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'results');

kernel_function = path(9:end);


end