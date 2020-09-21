function [setting] = getSetting(kernel_function,dist_type)

setting = struct('kernel',cell(1,1),'dist_type',cell(1,1),'nIters',cell(1,1),'st',cell(1,1),'ed',cell(1,1),'interval',cell(1,1),'r',cell(1,1));

path = 'Defaults/setting.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.kernel = kernel_function;
setting.dist_type = dist_type;
setting.st = str2double(INFO{2}{strcmp(INFO{1},'start:')==1});
setting.ed = str2double(INFO{2}{strcmp(INFO{1},'end:')==1});
setting.interval = str2double(INFO{2}{strcmp(INFO{1},'interval:')==1});
setting.nIters = str2double(INFO{2}{strcmp(INFO{1},'nIters:')==1});

setting.r = zeros(2,1);
setting.r(1) = str2double(INFO{2}{strcmp(INFO{1},'r_d11B:')==1});
setting.r(2) = str2double(INFO{2}{strcmp(INFO{1},'r_d18O:')==1});


end