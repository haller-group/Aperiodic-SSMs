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

F_I_net = [];

fig = figure

% fig = figure
% fontsize(fig, 20, "points")
% set(gcf,'color','w');
% box on
% grid on ;
% axis([0,500,-4,7])
% ax1 = gca;
% ax2 = axes( 'Position' ,[.2 .7 .35 .2]);
% box on
% grid on ;
NMTET = [];
% zo = [0.3*exp(-1i*pi/2),0.3*exp(1i*pi/2),0.5*exp(-1i*pi/2),0.5*exp(1i*pi/2),0.4*exp(1i*0.5),0.4*exp(-1i*0.5)];
zo = [6.2*exp(-1i*pi/2),1.2*exp(1i*pi/2),1*exp(-1i*pi/2),1*exp(1i*pi/2),0.8*exp(1i*0.5),0.1*exp(-1i*0.5)];
nssm = 30;
 x = linspace(-1,1,nssm);
     [X,Y]=meshgrid(x);
%      Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
      Net1 = zeros(nssm,nssm);
Net2 = zeros(nssm,nssm);

     Ratio_Net = zeros(nssm,nssm);
for ind = 1:nssm
    for indk = 1:nssm
epsilon = 0.01;
ctspan = linspace(0,6,4000)/epsilon;
ROM=@(t,z) rom_temp_model_adiabatic(t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_01,xi_11,xi_21,Valpha,V,A,Force_Lorenz,Dalpha,gamma,epsilon);
q0 = Z(ind,indk);
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,ctspan(1),gamma,q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

% tic 
% [t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
% toc
% 
% [y,modal,Net_Sol] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,t_sol,gamma,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);
% 

FullS = @(t,y) full_system(t,y,0,Force_Lorenz);
tic
IC =y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc
% % save('OrderOmega0.1.mat','SP_Traj','y','Net_Sol','tSP');
% 
yT = SP_Traj.';
% F_auto = A*[yT(1,:);yT(2,:);yT(3,:);yT(4,:);yT(5,:);yT(6,:)] + [0*yT(1,:);0*yT(1,:);0*yT(1,:);-gamma*(m1+Mf)/(m1*Mf)*yT(1,:).^3;gamma/m1*yT(1,:).^3;(Force_Lorenz(tSP*epsilon)/(m1+m2+Mf)).'];
% F_net = sqrt(sum(F_auto.^2));
% Ratio_Net(ind,indk) = epsilon./(trapz(tSP,F_net/max(F_net))/(max(tSP)));

epsilon = 0.008;
F_an = (abs(Force_Lorenz(epsilon*tSP))).';
% F_an = trapz(tSP,F_an)/(max(tSP)*mean(F_an));
F_spring = K*[yT(1,:);yT(2,:);yT(3,:)];
F_damper = C*[yT(4,:);yT(5,:);yT(6,:)];
F_nonlinear = [gamma*yT(1,:).^3;0*yT(1,:);0*yT(1,:)];
F_total = F_spring + F_damper + F_nonlinear;
F_net = sqrt(sum(F_total.^2));

[dXdF_net, Xtrunc_net,tSp_red] = ftd(F_net, tSP.');
[dXdF_an, Xtrunc_an,tSp_red] = ftd(F_an, tSP.');
F_I_net = [F_I_net;trapz(tSP,F_an./F_net)/(max(tSP))];

Net1(ind,indk) = trapz(tSp_red,abs(dXdF_an.'))/(max(tSP));
 Net2(ind,indk) = (trapz(tSp_red,abs(dXdF_net.'))/(max(tSP)));

Ratio_Net(ind,indk) = trapz(tSp_red,abs(dXdF_an.'))/(max(tSP))./(trapz(tSp_red,abs(dXdF_net.'))/(max(tSP)));


% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);

%  hold(ax1, 'on' );
% indexR = 3;
% plot(ax1,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
% hold on 
% plot(ax1,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
% hold on 
% plot(ax1,tSP,SolutionA(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
% hold on 
% xlabel(ax1,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax1,'$x_{c} \,[$m$]$','Interpreter','latex');
% legend(ax1,'Full Order Model','Reduced Order Model $O(5)$','Anchor Trajectory','Interpreter','latex')
% title(ax1,'\epsilon = 0.5, reduced order model performance')
% hold on
% 
% grid on ;
% 
% % ax2 = axes( 'Position' ,[.2 .7 .5 .2])
% hold(ax2, 'on' );
% hold on 
% plot(ax2,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
% hold on 
% plot(ax2,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
% hold on 
% plot(ax2,tSP,SolutionA(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
% hold on 
% xlabel(ax2,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax2,'$x_{c} \,[$m$]$','Interpreter','latex');
% axis([0,80,-4,4])
% grid on ;
% 
% NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));
% NMTET = [NMTET,NMTE];
    end
end
% origUnits = fig.Units;
% fig.Units = fig.PaperUnits;
% fig.PaperSize = fig.Position(3:4);
% fig.Units = origUnits;
% exportgraphics(fig, 'e008xcdot_slow.pdf');

% NMTEavg = sum(NMTET)/6;
fig = figure 
% flag = [Ratio_Net > 100];
% Ratio_Net(flag)= NaN;
bp = mean( Net1(:))./mean( Net2(:));
P = pcolor(X,Y,Ratio_Net);

% bp = mean(Ratio_Net(:));
P.EdgeColor = 'none';
colormap copper
colorbar 
xlabel('Re[$u_1$]','Interpreter','latex');
ylabel('Im[$u_1$]','Interpreter','latex');
title_string = strcat('$\displaystyle \frac{\overline{|\frac{\partial}{\partial t}|F_{ext}|_{t_0}^{t_f}(x_0)|}}{\overline{|\frac{d}{dt}|F_{int}|_{t_0}^{t_f}(x_0)|}}$, $\epsilon = ',num2str(epsilon),'$');


title(title_string,'Interpreter','latex')
daspect([1 1 1])

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'slow_small.pdf' , '-dpdf' , '-r300' )