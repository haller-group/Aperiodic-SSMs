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
a = 0.3;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
scale = 5*a^3;%2*a^3;
Sc = scale;
K = [-4*a^2*m*g/scale+k*m^2/(m+M)^2, -k*m/(m+M);-k*m/(m+M),k];
C = [cf+c*m^2/(m+M)^2,-c*m/(m+M);-c*m/(m+M),c];
MT = [m*M/(m+M),0;0,m+M];
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
% zo = [0.3*exp(-1i*pi/2),0.3*exp(1i*pi/2),0.5*exp(-1i*pi/2),0.5*exp(1i*pi/2),0.4*exp(1i*0.5),0.4*exp(-1i*0.5)];
% zo = [6.2*exp(-1i*pi/2),1.2*exp(1i*pi/2),1*exp(-1i*pi/2),1*exp(1i*pi/2),0.8*exp(1i*0.5),0.1*exp(-1i*0.5)];
nssm = 30;
 x = linspace(-1,1,nssm);
     [X,Y]=meshgrid(x);
%      Z = RHO.*exp(1i*Theta);
%      Z = X+1i*Y;
     Ratio_Net = zeros(nssm,nssm);
for ind = 1:nssm
    for indk = 1:nssm
epsilon = 1;

% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);
 ctspan = linspace(0,500,8000);
u0 = [X(ind,indk);Y(ind,indk)];
[y0,modal0] = compute_SSM_phy(XI_hyperbolic_solution,0,V_new,VI_new,u0,h030,h300,h120,h210,h111,h201,h021,h012,h102,epsilon);
GH = @(t,y) Full_Order_Model_PM_vec(t,y,Force_Lorenz,0,c,k,a,g,m,M,cf,Sc);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc
% q0 = Z(ind,indk);
% y0 =compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);
% 
% IC = [real(q0);imag(q0)];

% IC = [1.2;0];
% tic 
% [t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
% toc
% 
% lf = 0;
% [y,modal] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);




% SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
% FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);
% 
% tic
% IC =  y0;
% [tSP,SP_Traj] = ode45(FullS,ctspan,IC);
% toc

y = Full_Traj.';
F_an = (abs((m+M)*epsilon*Force_Lorenz(tSP))).';
% F_an = trapz(tSP,F_an)/(max(tSP)*mean(F_an));

F_spring = K*[y(1,:);y(2,:)];
F_damper = C*[y(3,:);y(4,:)];
F_nonlinear = [4*g*m*y(1,:).^3/scale + 16*m*a^4*y(1,:).*y(3,:).^2/scale^2;0*y(1,:)];
F_total = F_spring + F_damper + F_nonlinear;
F_net = sqrt(sum(F_total.^2));
F_I_net = [F_I_net;max(F_net)];

Ratio_Net(ind,indk) = trapz(tSP,F_an)/(max(tSP))./(trapz(tSP,F_net)/(max(tSP)));

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
factor = epsilon/max(F_I_net);
fig = figure
P = pcolor(X,Y,Ratio_Net);
P.EdgeColor = 'none';
colormap copper
colorbar 
xlabel('$u_1$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

title_string = strcat('$\displaystyle \frac{\overline{|F_{ext}|_{t_0}^{t_f}(x_0)}}{\overline{|F_{int}|_{t_0}^{t_f}(x_0)}}$, max$(|F_{ext}|) = ',num2str(epsilon*(M+m)),'$ [$\mathbf{N}$]');


title(title_string,'Interpreter','latex')

figure 
contour(X,Y,Ratio_Net)
colorbar