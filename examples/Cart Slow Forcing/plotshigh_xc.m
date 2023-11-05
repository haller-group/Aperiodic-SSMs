m1 = 1;
m2 = 1;
Mf = 2;
kf = 1;
k = 1;
c = 0.3;
cf = 0.3;
gamma = 0.5;
MTT = m1 + m2 + Mf;
K = [-2*k - kf*((m1 + m2)/MTT)^2, -k - kf*m2*(m1 + m2)/MTT^2, kf*(m1+m2)/MTT; -k - kf*m2*(m1 + m2)/MTT^2, -2*k - kf*m2^2/(MTT^2), kf*m2/MTT; (m1 + m2)/MTT *kf, m2/MTT*kf, -kf];
C = [-c - cf*((m1 + m2)/MTT)^2, -c - cf*m2*(m1 + m2)/MTT^2, cf*(m1+m2)/MTT; -c - cf*m2*(m1 + m2)/MTT^2, -2*c - cf*m2^2/(MTT^2), cf*m2/MTT; (m1 + m2)/MTT *cf, m2/MTT*cf, -cf];
M =  [Mf*(m1 + m2)/MTT, m2*Mf/MTT, 0; m2*Mf/MTT, m2*(m1 + Mf)/MTT, 0; 0,0, MTT];


A = [zeros(3),eye(3);inv(M)*K,inv(M)*C];
F_net = [];


fig = figure
fontsize(fig, 20, "points")
set(gcf,'color','w');
box on
grid on ;
axis([0,400,-11,18])
ax1 = gca;
ax2 = axes( 'Position' ,[.18 .7 .3 .2]);
box on
grid on ;
NMTET = [];
zo = [1.5*exp(-1i*pi/2),1.5*exp(1i*pi/2),2.5*exp(-1i*pi/2),2.5*exp(1i*pi/2),1*exp(1i*0.5),1*exp(-1i*0.5)];
for ind = 1:6
epsilon = 0.015;
ctspan = linspace(0,6,4000)/epsilon;
ROM=@(t,z) rom_temp_model_adiabatic(t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_01,xi_11,xi_21,Valpha,V,A,Force_Lorenz,Dalpha,gamma,epsilon);
q0 = zo(ind);
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,ctspan(1),gamma,q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

tic 
[t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
toc

[y,modal,Net_Sol] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,t_sol,gamma,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz);
tic
IC =y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc
% % save('OrderOmega0.1.mat','SP_Traj','y','Net_Sol','tSP');
% 
yT = SP_Traj.';
F_auto = A*[yT(1,:);yT(2,:);yT(3,:);yT(4,:);yT(5,:);yT(6,:)] + [0*yT(1,:);0*yT(1,:);0*yT(1,:);-gamma*(m1+Mf)/(m1*Mf)*yT(1,:).^3;gamma/m1*yT(1,:).^3;(Force_Lorenz(tSP*epsilon)/(m1+m2+Mf)).'];
F_net_1 = sqrt(sum(F_auto.^2));
F_net = [F_net;F_net_1];
 hold(ax1, 'on' );
indexR = 3;
plot(ax1,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ax1,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
hold on 
plot(ax1,tSP,Net_Sol(indexR,:),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
xlabel(ax1,'$t \,[$s$]$','Interpreter','latex');
ylabel(ax1,'$x_c \,[$m$]$','Interpreter','latex');
legend(ax1,'Full order model','Reduced order model','Chaotic anchor trajectory')
title(ax1,'$\epsilon = 0.008$, expansion order $N=3$','Interpreter','latex')
hold on

grid on ;

% ax2 = axes( 'Position' ,[.2 .7 .5 .2])
hold(ax2, 'on' );
hold on 
plot(ax2,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ax2,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
hold on 
plot(ax2,tSP,Net_Sol(indexR,:),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
% xlabel(ax2,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax2,'$x_{c} \,[$m$]$','Interpreter','latex');
axis([0,50,-10.5,1.5])
grid on ;

NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));
NMTET = [NMTET,NMTE];
end
% origUnits = fig.Units;
% fig.Units = fig.PaperUnits;
% fig.PaperSize = fig.Position(3:4);
% fig.Units = origUnits;
% exportgraphics(fig, 'e008xcdot_slow.pdf');

NMTEavg = sum(NMTET)/6;

% fig1 = figure
% fontsize(fig1, 20, "points")
% set(gcf,'color','w');
% box on
% grid on ;
% ax1 = gca;
% for ind = 1:6
%   plot(ax1,tSP,F_net(ind,:),'-','LineWidth',3,'color',[0 0 0 0.5])
% hold on 
%   
% end    
% hold on 
% plot(ax1,tSP,epsilon*ones(1,max(size(tSP))),'-','LineWidth',3,'color',[1 0 0 0.3])
% 
% save('NMTE0.008_3',"NMTEavg")
% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'cart_xc_slow0115.pdf' , '-dpdf' , '-r300' )
% 
% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'forcing_weak.pdf' , '-dpdf' , '-r300' )