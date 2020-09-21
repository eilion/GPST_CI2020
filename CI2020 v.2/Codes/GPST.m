function [results] = GPST(kernel_function,dist_type)

disp('-----------------------------------------------------------------');
disp('## This is a Gaussian process state-space model algorithm for reconstructing atmospheric CO2.');

% Load setting:
setting = getSetting(kernel_function,dist_type);
disp(['#  The kernel is ',setting.kernel,'.']);
disp(['#  The Observation model is ',dist_type,'.']);
% Load data:
[data,model,stack,setting,timeline] = getData(setting);
disp('#  The initialization step is done.');
disp('-----------------------------------------------------------------');

% Learn variational parameters and kernel hyperparameters:
disp(['## Reconstruction steps will be iterated for ',num2str(setting.nIters),' times.']);
[data,timeline] = learnParam(data,model,setting,timeline);
disp('   Done.');

% Variational Inference:
disp('## Variational Inference...');
[data,stack,setting] = varInference(data,stack,setting);
disp('   Done.');

results = struct('data',cell(1,1),'model',cell(1,1),'stack',cell(1,1),'setting',cell(1,1),'timeline',cell(1,1));
results.data = data;
results.model = model;
results.stack = stack;
results.setting = setting;
results.timeline = timeline;

disp('-----------------------------------------------------------------');
kernel_function = saveResults(kernel_function,dist_type,results);

getFigure(kernel_function);


end