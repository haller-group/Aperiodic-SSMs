clear all
m1 = 1;
m2 = 1;
Mf = 4;
kf = 1;
k = 1;
c = 0.3;
cf = 0.3;
gamma = 0.5;
MTT = m1 + m2 + Mf;
K = -[-2*k - kf*((m1 + m2)/MTT)^2, -k - kf*m2*(m1 + m2)/MTT^2, kf*(m1+m2)/MTT; -k - kf*m2*(m1 + m2)/MTT^2, -2*k - kf*m2^2/(MTT^2), kf*m2/MTT; (m1 + m2)/MTT *kf, m2/MTT*kf, -kf];
C = -[-c - cf*((m1 + m2)/MTT)^2, -c - cf*m2*(m1 + m2)/MTT^2, cf*(m1+m2)/MTT; -c - cf*m2*(m1 + m2)/MTT^2, -2*c - cf*m2^2/(MTT^2), cf*m2/MTT; (m1 + m2)/MTT *cf, m2/MTT*cf, -cf];
M =  [Mf*(m1 + m2)/MTT, m2*Mf/MTT, 0; m2*Mf/MTT, m2*(m1 + Mf)/MTT, 0; 0,0, MTT];

subs2 = [1 1 1 1];
n = 3;
F3 = sptensor(subs2, gamma, [n,n,n,n]);
F2 = sptensor([n,n,n]);
fnl = {F2,F3};

F_I_net = [];

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',6,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
S = SSM(DS);
resonant_modes = [1 2];
Acheck = V*diag(D)*inv(V);
VI = inv(V);
l1 = D(1);
l1c = conj(D(1));
l2 = D(3);
l2c = conj(D(3));
l3 = D(5);
l3c = conj(D(5));


[h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ] = Auto_Coeffs(gamma,l2,l3,l1c,l1);
load('Coeff_SSM_Exact_O1.mat','H_Coeff1','F_Coeff1');
load('Coeff_SSM_Exact_O2.mat','H_Coeff2','F_Coeff2');
load('Xi_Sol1.mat','xi_sol1')
load('Solution_3.mat','xi_full_eval3')
load('Solution_1.mat','xi_full_eval')
load('Solution_5.mat','xi_full_eval5')
load('Xi_Sol3.mat','xi_sol3')


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
 x = linspace(-2,2,nssm);
     [X,Y]=meshgrid(x);
%      Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
     Ratio_Net = zeros(nssm,nssm);
     accel = [];
for ind = 1:nssm
    for indk = 1:nssm
epsilon = 0.5;

% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);
 ctspan = linspace(0,500,8000);

q0 = Z(ind,indk);
y0 =compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);

IC = [real(q0);imag(q0)];

% IC = [1.2;0];
% tic 
% [t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
% toc
% 
% lf = 0;
% [y,modal] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);




% SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
FullS = @(t,y) full_system(t,y,0,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);

tic
IC =  y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc

y = SP_Traj.';
F_an = (abs((m1+m2+Mf)*epsilon*Force_Lorenz(tSP))).';
% F_an = trapz(tSP,F_an)/(max(tSP)*mean(F_an));

F_spring = K*[y(1,:);y(2,:);y(3,:)];
F_damper = C*[y(4,:);y(5,:);y(6,:)];
F_nonlinear = [gamma*y(1,:).^3;0*y(1,:);0*y(1,:)];
F_total = F_spring + F_damper + F_nonlinear;
F_net = sqrt(sum(F_total.^2));
F_I_net = [F_I_net;trapz(tSP,F_an./F_net)/(max(tSP))];
Ratio_Net(ind,indk) = trapz(tSP,F_an)/(max(tSP))./(trapz(tSP,F_net)/(max(tSP)));

F_spring = inv(M)*K*[y(1,:);y(2,:);y(3,:)];
F_damper = inv(M)*C*[y(4,:);y(5,:);y(6,:)];
F_nonlinear = inv(M)*[gamma*y(1,:).^3;0*y(1,:);0*y(1,:)];
F_total = F_spring + F_damper + F_nonlinear;
F_net = sqrt(sum(F_total.^2));
accel = [accel,F_net];

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
flag = [abs(sqrt(X.^2 + Y.^2)) < 0.3];
Ratio_Net(flag)= NaN;
NMTEavg = sum(NMTET)/6;
fig = figure
P = pcolor(X,Y,Ratio_Net);
P.EdgeColor = 'none';
colormap copper
colorbar 
xlabel('Re[$u_1$]','Interpreter','latex');
ylabel('Im[$u_1$]','Interpreter','latex');
title_string = strcat('$\displaystyle \frac{\overline{|F_{ext}|_{t_0}^{t_f}(x_0)}}{\overline{|F_{int}|_{t_0}^{t_f}(x_0)}}$, max$(|F_{ext}|) = ',num2str((m1+m2+Mf)*epsilon),'$ [$\mathbf{N}$]');
daspect([1 1 1])

title(title_string,'Interpreter','latex')
impratio = epsilon/max(accel);

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'cart_force_quantification_small.pdf' , '-dpdf' , '-r300' )
