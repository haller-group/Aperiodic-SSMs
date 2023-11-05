Epsilon = [0.1,0.05,0.01,0.007,0.005,0.003,0.001];
NMTET = [];
load('OrderOmega0.1.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p1 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.05.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p05 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.01.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p01 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.007.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p007 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.005.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p005 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.003.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p003 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

load('OrderOmega0.001.mat','SP_Traj','y','Net_Sol','tSP');
NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)));
NMTET = [NMTET,NMTE];
RelError0p001 = (sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

figure 
plot(Epsilon,NMTET,'-o','MarkerSize',30,'LineWidth',3,'color','blue')
xlabel('$\epsilon \,[$s^{-1}$]$','Interpreter','latex');
ylabel('$NMTE$','Interpreter','latex');

figure 
plot(tSP*0.001,RelError0p001,'-','LineWidth',3)
hold on
% plot(tSP*0.001,RelError0p003,'-','LineWidth',3)
% hold on 
% plot(tSP*0.001,RelError0p005,'-','LineWidth',3)
% hold on 
% plot(tSP*0.001,RelError0p007,'-','LineWidth',3)
% hold on 
% plot(tSP*0.001,RelError0p01,'-','LineWidth',3)
% hold on 
% plot(tSP*0.001,RelError0p05,'-','LineWidth',3)
% hold on 
% plot(tSP*0.001,RelError0p1,'-','LineWidth',3)

xlabel('$\alpha$','Interpreter','latex');
ylabel('Point-Wise Relative Error','Interpreter','latex');

