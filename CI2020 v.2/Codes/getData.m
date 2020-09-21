function [data,model,stack,setting,timeline] = getData(setting)

list = dir('Data/d11B/');
list = list(3:end);

L = length(list);

model = struct('d18O',cell(1,1),'d11B',cell(1,1));
model.d18O = struct('A',cell(1,1),'sig',cell(1,1));
model.d11B = struct('name',cell(L,1),'A',cell(L,1),'sig',cell(L,1));

for ll = 1:L
    model.d11B(ll).name = list(ll).name;
    
    path = ['Data/d11B/',list(ll).name,'/model.txt'];
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    model.d11B(ll).A = zeros(4,1);
    model.d11B(ll).A(1) = str2double(INFO{2}{strcmp(INFO{1},'a0:')==1});
    model.d11B(ll).A(2) = str2double(INFO{2}{strcmp(INFO{1},'a1:')==1});
    model.d11B(ll).A(3) = str2double(INFO{2}{strcmp(INFO{1},'a2:')==1});
    model.d11B(ll).A(4) = str2double(INFO{2}{strcmp(INFO{1},'a3:')==1});
    
    model.d11B(ll).sig = str2double(INFO{2}{strcmp(INFO{1},'sig:')==1});
end

path = 'Data/d18O_plankton/model.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

model.d18O.A = zeros(3,1);
model.d18O.A(1) = str2double(INFO{2}{strcmp(INFO{1},'a0:')==1});
model.d18O.A(2) = str2double(INFO{2}{strcmp(INFO{1},'a1:')==1});
model.d18O.A(3) = str2double(INFO{2}{strcmp(INFO{1},'a2:')==1});

model.d18O.sig = str2double(INFO{2}{strcmp(INFO{1},'sig:')==1});


data = struct('T1',cell(1,1),'T2',cell(1,1),'Y1',cell(1,1),'Y2',cell(1,1),'Z',cell(1,1),'mu1',cell(1,1),'sig1',cell(1,1),'mu2',cell(1,1),'sig2',cell(1,1),'eta',cell(1,1),'xi',cell(1,1),'lambda',cell(1,1),'delta',cell(1,1),'rho',cell(1,1),'rho0',cell(1,1));

data.eta = ones(2,1);
data.xi = 2*ones(3,1);
data.lambda = 0.1*ones(2,1);
data.lambda0 = 1e-4;

data.delta = ones(2,1);
data.rho0 = 0.3;
data.rho = tanh(data.rho0);

TT = cell(L,1);
YY = cell(L,1);
ZZ = cell(L,1);
for ll = 1:L
    path = ['Data/d11B/',model.d11B(ll).name,'/age.txt'];
    TT{ll} = load(path);
    
    path = ['Data/d11B/',model.d11B(ll).name,'/d11B.txt'];
    YY{ll} = load(path);
    
    N = size(TT{ll},1);
    ZZ{ll} = ll*ones(N,1);
end

data.T1 = cat(1,TT{:});
data.Y1 = cat(1,YY{:});
data.Z = cat(1,ZZ{:});

[~,order] = sort(data.T1,'ascend');
data.T1 = data.T1(order,:);
data.Y1 = data.Y1(order,:);
data.Z = data.Z(order,:);

stack = load('Data/d18O_plankton/stack.txt');

list = dir('Data/d18O_plankton/');
L = length(list);
index = 1:L;
for ll = 1:L
    if strcmp(list(ll).name,'.') == 1 || strcmp(list(ll).name,'..') == 1 || list(ll).isdir == 0
        index(ll) = NaN;
    end
end
index = index(~isnan(index));
list = list(index);

L = length(list);
TT = cell(L,1);
YY = cell(L,1);
for ll = 1:L
    path = ['Data/d18O_plankton/',list(ll).name,'/age.txt'];
    TT{ll} = load(path);
    
    path = ['Data/d18O_plankton/',list(ll).name,'/d18O.txt'];
    YY{ll} = load(path);
    
    h = mean(YY{ll}-interp1(stack(:,1),stack(:,2),TT{ll}));
    YY{ll} = YY{ll} - h;
end

data.T2 = cat(1,TT{:});
data.Y2 = cat(1,YY{:});

[~,order] = sort(data.T2,'ascend');
data.T2 = data.T2(order,:);
data.Y2 = data.Y2(order,:);

N1 = size(data.T1,1);
N2 = size(data.T2,1);

data.mu1 = 1e-2*normrnd(0,1,[N1,1]);
data.sig1 = 1e-4*(rand(N1,1)+1);

data.mu2 = 1e-2*normrnd(0,1,[N2,1]);
data.sig2 = 1e-4*(rand(N2,1)+1);

data.Y1 = (data.Y1-20.50)/1.5;


stack = struct('age',cell(1,1),'mu',cell(1,1),'sig',cell(1,1),'rho',cell(1,1));
stack.age = (setting.st:setting.interval:setting.ed)';

setting.mu_T = mean([data.T1;data.T2]);
setting.sig_T = sqrt(var([data.T1;data.T2]));

data.T1 = (data.T1-setting.mu_T)/setting.sig_T;
data.T2 = (data.T2-setting.mu_T)/setting.sig_T;
stack.age = (stack.age-setting.mu_T)/setting.sig_T;

setting.r = setting.r/setting.sig_T;


timeline = struct('delta',cell(1,1),'rho',cell(1,1),'eta',cell(1,1),'xi',cell(1,1),'lambda',cell(1,1),'abs_mean_mu',cell(1,1),'abs_mean_sig',cell(1,1));
timeline.delta = zeros(2,ceil(setting.nIters/100)+1);
timeline.rho = zeros(1,ceil(setting.nIters/100)+1);
timeline.eta = zeros(2,ceil(setting.nIters/100)+1);
timeline.xi = zeros(3,ceil(setting.nIters/100)+1);
timeline.lambda = zeros(2,ceil(setting.nIters/100)+1);
timeline.abs_mean_mu = zeros(1,ceil(setting.nIters/100)+1);
timeline.abs_mean_sig = zeros(1,ceil(setting.nIters/100)+1);

timeline.delta(:,1) = data.delta;
timeline.rho(1) = data.rho;
timeline.eta(:,1) = data.eta;
timeline.xi(:,1) = data.xi;
timeline.lambda(:,1) = data.lambda;
timeline.abs_mean_mu(1) = mean(abs([data.mu1;data.mu2]));
timeline.abs_mean_sig(1) = mean(abs([data.sig1;data.sig2]));


end