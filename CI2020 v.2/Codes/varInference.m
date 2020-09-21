function [data,stack,setting] = varInference(data,stack,setting)

N = size(stack.age,1);

N1 = size(data.T1,1);
N2 = size(data.T2,1);

MU = [data.mu1;data.mu2];
SIG = [data.sig1.^2;data.sig2.^2];

K_0_11 = getCov(data.T1,data.T1,data,setting,0);
K_0_12 = getCov(data.T1,data.T2,data,setting,0);
K_0_22 = getCov(data.T2,data.T2,data,setting,0);
K_1_11 = getCov(data.T1,data.T1,data,setting,1) + data.lambda(1)^2*eye(N1) + data.lambda0^2*eye(N1);
K_2_22 = getCov(data.T2,data.T2,data,setting,2) + data.lambda(2)^2*eye(N2) + data.lambda0^2*eye(N2);

KK = [data.delta(1)^2*K_0_11+K_1_11,data.rho*data.delta(1)*data.delta(2)*K_0_12;data.rho*data.delta(1)*data.delta(2)*K_0_12',data.delta(2)^2*K_0_22+K_2_22];

KK_1 = [data.delta(1)^2*getCov(stack.age,data.T1,data,setting,0)+getCov(stack.age,data.T1,data,setting,1),data.rho*data.delta(1)*data.delta(2)*getCov(stack.age,data.T2,data,setting,0)];
KK_2 = [data.rho*data.delta(1)*data.delta(2)*getCov(stack.age,data.T1,data,setting,0),data.delta(2)^2*getCov(stack.age,data.T2,data,setting,0)+getCov(stack.age,data.T2,data,setting,2)];

stack.mu = zeros(N,2);
stack.mu(:,1) = KK_1*(KK\MU);
stack.mu(:,2) = KK_2*(KK\MU);

EE1 = (KK\KK_1')';
EE2 = (KK\KK_2')';

stack.sig = zeros(N,2);
stack.sig(:,1) = sqrt(data.delta(1)^2 + data.eta(1)^2 + data.lambda(1)^2 - sum(EE1.*KK_1,2) + sum(EE1.*(EE1.*SIG'),2));
stack.sig(:,2) = sqrt(data.delta(2)^2 + data.eta(2)^2 + data.lambda(2)^2 - sum(EE2.*KK_2,2) + sum(EE2.*(EE2.*SIG'),2));

QQ = data.rho*data.delta(1)*data.delta(2) - sum(EE1.*KK_2,2) + sum(EE1.*(EE2.*SIG'),2);
stack.rho = QQ./(stack.sig(:,1).*stack.sig(:,2));


data.T1 = data.T1*setting.sig_T + setting.mu_T;
data.T2 = data.T2*setting.sig_T + setting.mu_T;
stack.age = stack.age*setting.sig_T + setting.mu_T;

data.Y1 = data.Y1*1.5 + 20.50;

setting.r = setting.r*setting.sig_T;

stack.mu(:,1) = stack.mu(:,1)*150 + 300;
stack.sig(:,1) = stack.sig(:,1)*150;


end