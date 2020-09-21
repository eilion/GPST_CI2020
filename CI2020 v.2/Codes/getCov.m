function [K] = getCov(X1,X2,data,setting,ll)

DIST_X = abs(X1(:,1)-X2(:,1)');

if ll == 0
    if strcmp(setting.kernel,'OU') == 1
        K = exp(-data.xi(1)^2*DIST_X);
    elseif strcmp(setting.kernel,'SE') == 1
        K = exp(-0.5*data.xi(1)^4*DIST_X.^2);
    elseif strcmp(setting.kernel,'M15') == 1
        K = (1+sqrt(3)*data.xi(1)^2*DIST_X).*exp(-sqrt(3)*data.xi(1)^2*DIST_X);
    elseif strcmp(setting.kernel,'M25') == 1
        K = (1+sqrt(5)*data.xi(1)^2*DIST_X+(sqrt(5)*data.xi(1)^2*DIST_X/sqrt(3)).^2).*exp(-sqrt(5)*data.xi(1)^2*DIST_X);
    end
else
    K = data.eta(ll)^2*exp(-2*data.xi(ll+1)^2*(sin(pi*DIST_X/setting.r(ll))).^2);
end


end