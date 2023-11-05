% m1 = 1;
% m2 = 1;
% Mf = 4;
% kf = 1;
% k = 1;
% c = 0.3;
% cf = 0.3;
% gamma = 0.5;
% MTT = m1 + m2 + Mf;
% K = -[-2*k - kf*((m1 + m2)/MTT)^2, -k - kf*m2*(m1 + m2)/MTT^2, kf*(m1+m2)/MTT; -k - kf*m2*(m1 + m2)/MTT^2, -2*k - kf*m2^2/(MTT^2), kf*m2/MTT; (m1 + m2)/MTT *kf, m2/MTT*kf, -kf];
% C = -[-c - cf*((m1 + m2)/MTT)^2, -c - cf*m2*(m1 + m2)/MTT^2, cf*(m1+m2)/MTT; -c - cf*m2*(m1 + m2)/MTT^2, -2*c - cf*m2^2/(MTT^2), cf*m2/MTT; (m1 + m2)/MTT *cf, m2/MTT*cf, -cf];
% M =  [Mf*(m1 + m2)/MTT, m2*Mf/MTT, 0; m2*Mf/MTT, m2*(m1 + Mf)/MTT, 0; 0,0, MTT];
% F_I_net = [];

c = 0.3;
k = 1;
a = 0.5;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
scale = 10*a^3;%2*a^3;
Sc = scale;
K = -[-k,k*m/(m+Mf);k*m/(m+Mf),4*a^2*m*g/scale-k*m^2/(m+Mf)^2];
C = -[-c,c*m/(m+Mf);c*m/(m+Mf),-cf-c*m^2/(m+Mf)^2];
M = [m+Mf,0;0,m*Mf/(m+Mf)];

F_I_net = [];

fig = figure
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
load('Solution_Unstable.mat','XI_0U','XI_1U','XI_2U','XI_3U');
load('Solution_stablem.mat','XI_0Sm','XI_1Sm','XI_2Sm','XI_3Sm');
load('Solution_stablep.mat','XI_0Sp','XI_1Sp','XI_2Sp','XI_3Sp');


% zo = [0.3*exp(-1i*pi/2),0.3*exp(1i*pi/2),0.5*exp(-1i*pi/2),0.5*exp(1i*pi/2),0.4*exp(1i*0.5),0.4*exp(-1i*0.5)];
% zo = [6.2*exp(-1i*pi/2),1.2*exp(1i*pi/2),1*exp(-1i*pi/2),1*exp(1i*pi/2),0.8*exp(1i*0.5),0.1*exp(-1i*0.5)];
nssm = 20;
 x = linspace(-1,1,nssm);
     [X,Y]=meshgrid(x);
%      Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
     Ratio_Net = zeros(nssm,nssm);
for ind = 1:nssm
    for indk = 1:nssm

epsilon = 0.01;

alphaT = linspace(0,6,10000);



ctspan = alphaT/(epsilon);



ROM=@(t,z) rom_temp_model_adiabatic_rail(t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_20,xi_21,xi_41,xi_22,xi_42,Valpha,V,A,Force_Lorenz,Dalpha,epsilon,a,g,c,cf,Mf,m,k,Sc);

q0 = Z(ind,indk);
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,ctspan(1),q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

% IC = [real(q0);imag(q0)];
% 
% tic 
% [t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
% toc

% [y,modal,Net_Sol] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_sol,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,0,c,k,a,g,m,Mf,cf,scale);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc







% SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
% FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);
% 
% tic
% IC =  y0;
% [tSP,SP_Traj] = ode45(FullS,ctspan,IC);
% toc

y = Full_Traj.';
F_an = (abs(Force_Lorenz(epsilon*tSP))).';
% F_an = trapz(tSP,F_an)/(max(tSP)*mean(F_an));

F_spring = K*[y(1,:);y(2,:)];
F_damper = C*[y(3,:);y(4,:)];
F_nonlinear = [0*y(1,:);4*g*m*y(2,:).^3/scale + 16*m*a^4*y(2,:).*y(4,:).^2/scale^2];
F_total = F_spring + F_damper + F_nonlinear;
F_net = sqrt(sum(F_total.^2));
[dXdF_net, Xtrunc_net,tSp_red] = ftd(F_net, tSP.');
[dXdF_an, Xtrunc_an,tSp_red] = ftd(F_an, tSP.');
F_I_net = [F_I_net;trapz(tSP,F_an./F_net)/(max(tSP))];

Ratio_Net(ind,indk) = trapz(tSp_red,abs(dXdF_an.'))/(max(tSP))./(trapz(tSp_red,abs(dXdF_net.'))/(max(tSP)));

% Ratio_Net(ind,indk) = trapz(tSP,F_an./F_net)/(max(tSP));

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

bp = mean(Ratio_Net(:));
NMTEavg = sum(NMTET)/6;

fig = figure
P = pcolor(X,Y,Ratio_Net);
P.EdgeColor = 'none';
colormap copper
colorbar 
xlabel('$u_1$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

title_string = strcat('$\displaystyle \frac{\overline{|\frac{d}{dt}|F_{ext}|_{t_0}^{t_f}(x_0)|}}{\overline{|\frac{d}{dt}|F_{int}|_{t_0}^{t_f}(x_0)|}}$, $\epsilon = ',num2str(epsilon),'$');


title(title_string,'Interpreter','latex')

figure 
contour(X,Y,Ratio_Net)
colorbar