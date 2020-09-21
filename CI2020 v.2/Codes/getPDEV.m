function [PDEV,PDEV1,PDEV2] = getPDEV(data,model,setting,rand_seed_1,rand_seed_2)

N1 = size(data.T1,1);
N2 = size(data.T2,1);

MU = [data.mu1;data.mu2];
SIG = [data.sig1.^2;data.sig2.^2];

PDEV = zeros(10,1);
PDEV1 = zeros(N1,2);
PDEV2 = zeros(N2,2);

K_0_11 = getCov(data.T1,data.T1,data,setting,0);
K_0_12 = getCov(data.T1,data.T2,data,setting,0);
K_0_22 = getCov(data.T2,data.T2,data,setting,0);
K_1_11 = getCov(data.T1,data.T1,data,setting,1) + data.lambda(1)^2*eye(N1) + data.lambda0^2*eye(N1);
K_2_22 = getCov(data.T2,data.T2,data,setting,2) + data.lambda(2)^2*eye(N2) + data.lambda0^2*eye(N2);

KK = [data.delta(1)^2*K_0_11+K_1_11,data.rho*data.delta(1)*data.delta(2)*K_0_12;data.rho*data.delta(1)*data.delta(2)*K_0_12',data.delta(2)^2*K_0_22+K_2_22];
KK_inv = KK\eye(N1+N2);
QQ = KK_inv*(diag(SIG)+MU*MU'-KK)*KK_inv;

% Kernel hyperparameters:
% delta:
PK = [2*data.delta(1)*K_0_11,data.rho*data.delta(2)*K_0_12;data.rho*data.delta(2)*K_0_12',zeros(N2)];
PDEV(8) = 0.5*sum(sum(QQ.*PK,2));

PK = [zeros(N1),data.rho*data.delta(1)*K_0_12;data.rho*data.delta(1)*K_0_12',2*data.delta(2)*K_0_22];
PDEV(9) = 0.5*sum(sum(QQ.*PK,2));

%{
% rho:
PK = [zeros(N1),data.delta(1)*data.delta(2)*K_0_12;data.delta(1)*data.delta(2)*K_0_12',zeros(N2)];
PDEV(10) = 0.5*sum(sum(QQ.*PK,2));
%}
% rho0:
PK = [zeros(N1),data.delta(1)*data.delta(2)*K_0_12;data.delta(1)*data.delta(2)*K_0_12',zeros(N2)]/(cosh(data.rho0)^2);
PDEV(10) = 0.5*sum(sum(QQ.*PK,2));

% eta:
PK_1_11 = 2*getCov(data.T1,data.T1,data,setting,1)/data.eta(1);
PK = [PK_1_11,zeros(N1,N2);zeros(N2,N1),zeros(N2)];
PDEV(1) = 0.5*sum(sum(QQ.*PK,2));

PK_2_22 = 2*getCov(data.T2,data.T2,data,setting,2)/data.eta(2);
PK = [zeros(N1),zeros(N1,N2);zeros(N2,N1),PK_2_22];
PDEV(2) = 0.5*sum(sum(QQ.*PK,2));

% xi:
DIST = abs(data.T1-data.T1');
if strcmp(setting.kernel,'OU') == 1
    PK_0_11 = exp(-data.xi(1)^2.*DIST).*(-2*data.xi(1).*DIST);
elseif strcmp(setting.kernel,'SE') == 1
    PK_0_11 = exp(-0.5*data.xi(1)^4.*DIST.^2).*(-2*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M15') == 1 || strcmp(setting.kernel,'M15PR') == 1
    PK_0_11 = exp(-sqrt(3)*data.xi(1)^2.*DIST).*(-2*sqrt(3)*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M25') == 1
    PK_0_11 = exp(-sqrt(5)*data.xi(1)^2.*DIST).*(-10*data.xi(1)^3.*DIST.^2)/3.*(1+sqrt(5)*data.xi(1)^2.*DIST);    
end
PK_1_11 = data.eta(1)^2*exp(-2*data.xi(2)^2*(sin(pi*DIST/setting.r(1))).^2).*(-4*data.xi(2)*(sin(pi*DIST/setting.r(1))).^2);

DIST = abs(data.T1-data.T2');
if strcmp(setting.kernel,'OU') == 1
    PK_0_12 = exp(-data.xi(1)^2.*DIST).*(-2*data.xi(1).*DIST);
elseif strcmp(setting.kernel,'SE') == 1
    PK_0_12 = exp(-0.5*data.xi(1)^4.*DIST.^2).*(-2*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M15') == 1 || strcmp(setting.kernel,'M15PR') == 1
    PK_0_12 = exp(-sqrt(3)*data.xi(1)^2.*DIST).*(-2*sqrt(3)*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M25') == 1
    PK_0_12 = exp(-sqrt(5)*data.xi(1)^2.*DIST).*(-10*data.xi(1)^3.*DIST.^2)/3.*(1+sqrt(5)*data.xi(1)^2.*DIST);    
end

DIST = abs(data.T2-data.T2');
if strcmp(setting.kernel,'OU') == 1
    PK_0_22 = exp(-data.xi(1)^2.*DIST).*(-2*data.xi(1).*DIST);
elseif strcmp(setting.kernel,'SE') == 1
    PK_0_22 = exp(-0.5*data.xi(1)^4.*DIST.^2).*(-2*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M15') == 1 || strcmp(setting.kernel,'M15PR') == 1
    PK_0_22 = exp(-sqrt(3)*data.xi(1)^2.*DIST).*(-2*sqrt(3)*data.xi(1)^3.*DIST.^2);
elseif strcmp(setting.kernel,'M25') == 1
    PK_0_22 = exp(-sqrt(5)*data.xi(1)^2.*DIST).*(-10*data.xi(1)^3.*DIST.^2)/3.*(1+sqrt(5)*data.xi(1)^2.*DIST);    
end
PK_2_22 = data.eta(2)^2*exp(-2*data.xi(3)^2*(sin(pi*DIST/setting.r(2))).^2).*(-4*data.xi(3)*(sin(pi*DIST/setting.r(2))).^2);

PK = [data.delta(1)^2*PK_0_11,data.rho*data.delta(1)*data.delta(2)*PK_0_12;data.rho*data.delta(1)*data.delta(2)*PK_0_12',data.delta(2)^2*PK_0_22];
PDEV(3) = 0.5*sum(sum(QQ.*PK,2));

PK = [PK_1_11,zeros(N1,N2);zeros(N2,N1),zeros(N2)];
PDEV(4) = 0.5*sum(sum(QQ.*PK,2));

PK = [zeros(N1),zeros(N1,N2);zeros(N2,N1),PK_2_22];
PDEV(5) = 0.5*sum(sum(QQ.*PK,2));

PK = [2*data.lambda(1)*eye(N1),zeros(N1,N2);zeros(N2,N1),zeros(N2)];
PDEV(6) = 0.5*sum(sum(QQ.*PK,2));

PK = [zeros(N1),zeros(N1,N2);zeros(N2,N1),2*data.lambda(2)*eye(N2)];
PDEV(7) = 0.5*sum(sum(QQ.*PK,2));


% Variational Parameters:
X1 = data.mu1 + data.sig1.*rand_seed_1;
X2 = data.mu2 + data.sig2.*rand_seed_2;

MU1 = zeros(N1,100);
for m = 1:length(model.d11B)
    index = (data.Z==m);
    MU1(index,:) = model.d11B(m).A(1) + model.d11B(m).A(2)*X1(index,:) + model.d11B(m).A(3)*log(model.d11B(m).A(4)+X1(index,:));
end
MU2 = model.d18O.A(1) + model.d18O.A(2)*X2 + model.d18O.A(3)*X2.^2;

QQ1 = zeros(N1,100);
for m = 1:length(model.d11B)
    index = (data.Z==m);
    if strcmp(setting.dist_type,'Normal') == 1
        QQ1(index,:) = model.d11B(m).sig^(-2)*(data.Y1(index)-MU1(index,:)).*(model.d11B(m).A(2)+model.d11B(m).A(3)./(model.d11B(m).A(4)+X1(index,:)));
    elseif strcmp(setting.dist_type,'T') == 1
        QQ1(index,:) = 7./(8+((data.Y1(index)-MU1(index,:))/model.d11B(m).sig).^2).*model.d11B(m).sig.^(-2).*(data.Y1(index)-MU1(index,:)).*(model.d11B(m).A(2)+model.d11B(m).A(3)./(model.d11B(m).A(4)+X1(index,:)));
    end
end
if strcmp(setting.dist_type,'Normal') == 1
    QQ2 = model.d18O.sig^(-2)*(data.Y2-MU2).*(model.d18O.A(2)+2*model.d18O.A(3)*X2);
elseif strcmp(setting.dist_type,'T') == 1
    QQ2 = 7./(8*model.d18O.sig.^2+(data.Y2-MU2).^2).*(data.Y2-MU2).*(model.d18O.A(2)+2*model.d18O.A(3)*X2);
end

% mu:
AA = KK_inv*MU;
PDEV1(:,1) = mean(QQ1,2) - AA(1:N1);
PDEV2(:,1) = mean(QQ2,2) - AA(N1+1:end);

% sig:
PSIG = [2*data.sig1;2*data.sig2];
AA = sum(diag(KK_inv).*PSIG,2);
PDEV1(:,2) = mean(QQ1.*rand_seed_1,2) - 0.5*AA(1:N1) + 1./data.sig1;
PDEV2(:,2) = mean(QQ2.*rand_seed_2,2) - 0.5*AA(N1+1:end) + 1./data.sig2;


end