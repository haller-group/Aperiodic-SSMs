
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
F_I_net = [];


fig = figure
fontsize(fig, 20, "points")
set(gcf,'color','w');
box on
grid on ;
axis([0,500,-4,7])
ax1 = gca;
ax2 = axes( 'Position' ,[.18 .7 .3 .2]);
box on
grid on ;
NMTET = [];
% zo = [0.3*exp(-1i*pi/2),0.3*exp(1i*pi/2),0.5*exp(-1i*pi/2),0.5*exp(1i*pi/2),0.4*exp(1i*0.5),0.4*exp(-1i*0.5)];
% zo = [1.2*exp(-1i*pi/2),1.2*exp(1i*pi/2),1*exp(-1i*pi/2),1*exp(1i*pi/2),0.8*exp(1i*0.5),0.8*exp(-1i*0.5)];
zo = [1.2*exp(-1i*pi/2),1.2*exp(1i*pi/2),0.8*exp(-1i*pi/2),0.8*exp(1i*pi/2),1*exp(1i*0.5),1*exp(-1i*0.5)];

for ind = 1:6
epsilon = 0.5;

ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);
ctspan = linspace(0,500,8000);

q0 = zo(ind);
y0 =compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);

IC = [real(q0);imag(q0)];

% IC = [1.2;0];
tic 
[t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
toc

lf = 0;
[y,modal] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);




SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
FullS = @(t,y) full_system(t,y,epsilon,Force_Dis_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);

tic
IC =  y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc

F_spring = K*[y(1,:);y(2,:);y(3,:)];
F_damper = C*[y(4,:);y(5,:);y(6,:)];
F_total = F_spring+F_damper;
F_net = sqrt(sum(F_total.^2));
F_I_net = [F_I_net;F_net];

 hold(ax1, 'on' );
indexR = 3;
plot(ax1,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ax1,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
hold on 
plot(ax1,tSP,SolutionA(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
xlabel(ax1,'$t \,[$s$]$','Interpreter','latex');
ylabel(ax1,'$x_{c} \,[$m$]$','Interpreter','latex');
legend(ax1,'Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex')
title(ax1,'max$(|F_{ext}|) = 3$ [$\mathbf{N}$], expansion order $N = 5$','Interpreter','latex')

hold on

grid on ;

% ax2 = axes( 'Position' ,[.2 .7 .5 .2])
hold(ax2, 'on' );
hold on 
plot(ax2,tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ax2,tSP,y(indexR,:),'--','LineWidth',3,'color','red')
hold on 
plot(ax2,tSP,SolutionA(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
% xlabel(ax2,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax2,'$x_{c} \,[$m$]$','Interpreter','latex');
xticks(ax2,[0 20 40])
yticks(ax2,[-3 0 2])
axis([0,40,-3.5,3])
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

% NMTEavg = 0.0063
% save('NMTE0.008_3',"NMTEavg")

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'cart_xc_weak_discontinous.pdf' , '-dpdf' , '-r300' )

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'discontinous_forcing.pdf' , '-dpdf' , '-r300' )
