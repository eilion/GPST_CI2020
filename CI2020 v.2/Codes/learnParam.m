function [data,timeline] = learnParam(data,model,setting,timeline)

N1 = length(data.T1);
N2 = length(data.T2);

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 2.5*1e-4;

max_iters = setting.nIters;

Mw = zeros(10,1);
Vw = zeros(10,1);

Mw1 = zeros(N1,2);
Vw1 = zeros(N1,2);

Mw2 = zeros(N2,2);
Vw2 = zeros(N2,2);

r = 0;
rand_seed_1 = normrnd(0,1,[N1,100]);
rand_seed_2 = normrnd(0,1,[N2,100]);
while r < max_iters
    r = r + 1;
    
    [PDEV,PDEV1,PDEV2] = getPDEV(data,model,setting,rand_seed_1,rand_seed_2);
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
    
    Mw1 = beta1*Mw1 - (1-beta1)*PDEV1;
    Vw1 = beta2*Vw1 + (1-beta2)*PDEV1.*PDEV1;
    
    Mw2 = beta1*Mw2 - (1-beta1)*PDEV2;
    Vw2 = beta2*Vw2 + (1-beta2)*PDEV2.*PDEV2;
    
    data.delta = abs(data.delta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(8:9)./(sqrt(Vw(8:9))+epsilon));
    data.rho0 = data.rho0 - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(10)./(sqrt(Vw(10))+epsilon);
    data.rho = tanh(data.rho0);
    
    data.eta = abs(data.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1:2)./(sqrt(Vw(1:2))+epsilon));
    data.xi = abs(data.xi - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3:5)./(sqrt(Vw(3:5))+epsilon));
    data.lambda = abs(data.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(6:7)./(sqrt(Vw(6:7))+epsilon));
    
    data.mu1 = data.mu1 - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw1(:,1)./(sqrt(Vw1(:,1))+epsilon);
    data.sig1 = abs(data.sig1 - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw1(:,2)./(sqrt(Vw1(:,2))+epsilon));
    
    data.mu2 = data.mu2 - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw2(:,1)./(sqrt(Vw2(:,1))+epsilon);
    data.sig2 = abs(data.sig2 - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw2(:,2)./(sqrt(Vw2(:,2))+epsilon));
    
    if rem(r,100) == 0
        timeline.delta(:,1+r/100) = data.delta;
        timeline.rho(1+r/100) = data.rho;
        timeline.eta(:,1+r/100) = data.eta;
        timeline.xi(:,1+r/100) = data.xi;
        timeline.lambda(:,1+r/100) = data.lambda;
        timeline.abs_mean_mu(1+r/100) = mean(abs([data.mu1;data.mu2]));
        timeline.abs_mean_sig(1+r/100) = mean(abs([data.sig1;data.sig2]));
    end
    
    if rem(r,1000) == 0
        disp(['#  Iterations ',num2str(r),'/',num2str(max_iters),' are done.']);
    end
end


end