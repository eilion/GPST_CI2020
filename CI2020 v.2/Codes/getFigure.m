function getFigure(filepath)

path = ['Outputs/',filepath,'/results.mat'];
results = load(path);
results = results.results;

stack = results.stack;
data = results.data;
setting = results.setting;
timeline = results.timeline;

bubble_CO2 = load('Benchmark/Bubble_CO2.txt');
published = load('Benchmark/Published.txt');
SST_stack = load('Benchmark/stack.txt');
MCMC = load('Benchmark/MCMC.txt');


fig = figure;

subplot(3,1,1);
h = zeros(3,1);
hold on;
if strcmp(setting.dist_type,'Normal')
    title(['Atmospheric CO2 (',setting.kernel,', Gaussian Observation)'],'FontSize',16);
elseif strcmp(setting.dist_type,'T')
    title(['Atmospheric CO2 (',setting.kernel,', T Observation)'],'FontSize',16);
end
for n = 1:size(published,1)
    h(2) = plot([published(n,1),published(n,1)],[published(n,2)+published(n,3),published(n,2)-published(n,4)],'c','LineWidth',2);
end
xx = [stack.age;flipud(stack.age)];
yy = [stack.mu(:,1)-2*stack.sig(:,1);flipud(stack.mu(:,1)+2*stack.sig(:,1))];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
plot(stack.age,stack.mu(:,1)+2*stack.sig(:,1),':k','LineWidth',2);
plot(stack.age,stack.mu(:,1)-2*stack.sig(:,1),':k','LineWidth',2);
h(3) = plot(bubble_CO2(:,1),bubble_CO2(:,2),'r','LineWidth',2);
h(1) = plot(stack.age,stack.mu(:,1),'--k','LineWidth',2);

xlim([setting.st setting.ed]);
ylim([100 450]);

xlabel('ages (kyr)','FontSize',12);
ylabel('CO_{2} (ppmv)','FontSize',12);

% legend(h,{'inferred','published','ice core'},'Location','NorthEast');

subplot(3,1,2);
h = zeros(3,1);
hold on;
title('Sea Surface Temperature Index','FontSize',16);

for n = 1:size(MCMC,1)
    h(2) = plot([MCMC(n,1),MCMC(n,1)],[MCMC(n,2),MCMC(n,3)],'c','LineWidth',2);
end
xx = [stack.age;flipud(stack.age)];
yy = [stack.mu(:,2)-2*stack.sig(:,2);flipud(stack.mu(:,2)+2*stack.sig(:,2))];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
plot(stack.age,stack.mu(:,2)+2*stack.sig(:,2),':k','LineWidth',2);
plot(stack.age,stack.mu(:,2)-2*stack.sig(:,2),':k','LineWidth',2);
h(3) = plot(SST_stack(:,1),SST_stack(:,2),'r','LineWidth',2);
h(1) = plot(stack.age,stack.mu(:,2),'--k','LineWidth',2);
% h(2) = plot(data.T2,-4.5*ones(size(data.T2,1),1),'^c','LineWidth',1);

xlim([setting.st setting.ed]);
ylim([-5 5]);

xlabel('ages (kyr)','FontSize',12);
ylabel('SST indices','FontSize',12);

% legend(h,{'inferred','individual','SST stack'},'Location','NorthEast');

subplot(3,1,3);
hold on;
title('Inferred Correlation','FontSize',16);
plot(stack.age,stack.rho,'k','LineWidth',2);
xlim([setting.st setting.ed]);
ylim([-1 1]);

xlabel('ages (kyr)','FontSize',12);
ylabel('correlation','FontSize',12);

% set(fig,'Position',[20 20 2000 1200]);
set(fig,'Position',[20 20 1000 900]);
movegui(fig,'center');
drawnow;

path = ['Outputs/',filepath,'/Results.fig'];
savefig(fig,path);



fig = figure;

subplot(2,3,1);
h = zeros(2,1);
hold on;
title('\delta','FontSize',16);
h(1) = plot(0:100:ceil(setting.nIters/100)*100,timeline.delta(1,:),'r','LineWidth',2);
h(2) = plot(0:100:ceil(setting.nIters/100)*100,timeline.delta(2,:),'g','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
ylim([0 max(max(timeline.delta))*1.05]);
legend(h,{'{\delta}_1','{\delta}_2'},'Location','East');

subplot(2,3,2);
hold on;
title('\rho','FontSize',16);
plot(0:100:ceil(setting.nIters/100)*100,timeline.rho,'k','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
% ylim([-max(abs(timeline.rho))*1.05 max(abs(timeline.rho))*1.05]);
ylim([-1 1]);

subplot(2,3,3);
h = zeros(2,1);
hold on;
title('\eta','FontSize',16);
h(1) = plot(0:100:ceil(setting.nIters/100)*100,timeline.eta(1,:),'r','LineWidth',2);
h(2) = plot(0:100:ceil(setting.nIters/100)*100,timeline.eta(2,:),'g','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
ylim([0 max(max(timeline.eta))*1.05]);
legend(h,{'{\eta}_1','{\eta}_2'},'Location','NorthEast');

subplot(2,3,4);
h = zeros(3,1);
hold on;
title('\xi','FontSize',16);
h(1) = plot(0:100:ceil(setting.nIters/100)*100,timeline.xi(1,:),'k','LineWidth',2);
h(2) = plot(0:100:ceil(setting.nIters/100)*100,timeline.xi(2,:),'r','LineWidth',2);
h(3) = plot(0:100:ceil(setting.nIters/100)*100,timeline.xi(3,:),'g','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
ylim([0 max(max(timeline.xi))*1.05]);
legend(h,{'{\xi}_0','{\xi}_1','{\xi}_2'},'Location','East');

subplot(2,3,5);
h = zeros(2,1);
hold on;
title('\lambda','FontSize',16);
h(1) = plot(0:100:ceil(setting.nIters/100)*100,timeline.lambda(1,:),'r','LineWidth',2);
h(2) = plot(0:100:ceil(setting.nIters/100)*100,timeline.lambda(2,:),'g','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
ylim([0 max(max(timeline.lambda))*1.05]);
legend(h,{'{\lambda}_1','{\lambda}_2'},'Location','East');

subplot(2,3,6);
h = zeros(2,1);
hold on;
title('Averages of Abs. Var. Param.','FontSize',16);
h(1) = plot(0:100:ceil(setting.nIters/100)*100,timeline.abs_mean_mu,'m','LineWidth',2);
h(2) = plot(0:100:ceil(setting.nIters/100)*100,timeline.abs_mean_sig,'c','LineWidth',2);
xlabel('iteration','FontSize',12);
xlim([0 ceil(setting.nIters/100)*100]);
ylim([0 max(max(timeline.abs_mean_mu),max(timeline.abs_mean_sig))*1.05]);
legend(h,{'mu','sig'},'Location','SouthEast');


set(fig,'Position',[20 20 1000 600]);
movegui(fig,'center');
drawnow;

path = ['Outputs/',filepath,'/Hyperparameters.fig'];
savefig(fig,path);


close all;


end