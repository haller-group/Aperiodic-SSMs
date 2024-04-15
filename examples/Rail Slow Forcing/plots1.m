

fig = figure
fontsize(fig, 20, "points")
set(gcf,'color','w');
box on
grid on ;
axis([0,400,-0.7,1.4])
ax1 = gca;
ax2 = axes( 'Position' ,[.2 .7 .35 .2]); 
box on
grid on ;
axis([0,600,0.46,0.54])
ax3 = axes( 'Position' ,[.65 .7 .2 .2]); 
box on
grid on ;
axis([150,152,-0.7,0.7])

NMTET = [];
% zo = [0.1*exp(1i*0),0.1*exp(1i*pi/6),0.1*exp(1i*pi/3),0.1*exp(1i*pi/2),0.1*exp(2*1i*pi/3),0.1*exp(5*1i*pi/6)];

b = 0.3;
zo = [ b+1i*0, ...
       -b+1i*0, ...
       b+1i*-b, ...
       b+1i*b,...
       -b+1i*-b,...
       -b+b*1i,...
       0+1i*b,...
       0+1i*-b,...
       2*b+1i*0, ...
       -2*b+1i*0,...
       2*b+1i*-2*b,...
       2*b+1i*2*b,...
       -2*b+1i*-2*b,...
       -2*b+2*b*1i,...
       0+1i*2*b,...
       0+1i*-2*b
    ];

epsilon = 0.01;

alphaTOG = linspace(0,6,100000);


to1=find(alphaTOG>=0.5);
to2=find(alphaTOG>=1);
to3=find(alphaTOG>=1.5);

to = [1,to1(1),to2(1),to3(1)];

XIU =[];
XISP=[];
XISM=[];

load('Solution_Unstable.mat','XI_0U','XI_1U','XI_2U','XI_3U');
load('Solution_stablem.mat','XI_0Sm','XI_1Sm','XI_2Sm','XI_3Sm');
load('Solution_stablep.mat','XI_0Sp','XI_1Sp','XI_2Sp','XI_3Sp');

for ind = 1:max(size(alphaTOG))
    alpha=alphaTOG(ind);
XI_Unstable = XI_0U(alpha)+ epsilon*XI_1U(alpha)+epsilon^2*XI_2U(alpha)+epsilon^3*XI_3U(alpha);
XIU =[XIU;XI_Unstable];
XI_Stablep = XI_0Sp(alpha)+ epsilon*XI_1Sp(alpha)+epsilon^2*XI_2Sp(alpha)+epsilon^3*XI_3Sp(alpha)+ [a*m/(m+Mf),a,0,0];
XISP =[XISP;XI_Stablep];
XI_Stablem = XI_0Sm(alpha)+ epsilon*XI_1Sm(alpha)+epsilon^2*XI_2Sm(alpha)+epsilon^3*XI_3Sm(alpha)- [a*m/(m+Mf),a,0,0];
XISM =[XISM;XI_Stablem];
end

store_Traj = cell(5,4);
store_Traj(1,:) = {alphaTOG/epsilon,XIU,XISP,XISM};
for indk = 1:4
Yrec = [];
FullTra = [];
    for ind = 1:16

epsilon = 0.01;

alphaT = alphaTOG(to(indk):end);
ctspan = alphaT/(epsilon);
stop_index = 500;

ctspanR = ctspan(1:stop_index);

ROM=@(t,z) rom_temp_model_adiabatic_rail(t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_20,xi_21,xi_41,xi_22,xi_42,Valpha,ValphaD,V,A,Force_Lorenz,Dalpha,epsilon,a,g,c,cf,Mf,m,k,Sc);

q0 = zo(ind);
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,ctspan(1),q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

tic 
[t_sol,ROM_sol] = ode45(ROM,ctspanR,IC);
toc

[y,modal,Net_Sol] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_sol,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon,c,k,a,g,m,Mf,cf,scale);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc


% [y0,model0,Net_Sol0] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,ctspan(1),q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);
% 
% IC = [real(q0);imag(q0)];
% 
% tic 
% [t_sol,ROM_sol] = ode45(ROM,ctspanR,IC);
% toc
% 
% [y,modal,Net_Sol] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_sol,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);
% 
% 
% GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon);
% tic
% IC = y0;
% [tSP,Full_Traj] = ode45(GH,ctspan,IC);
% toc




% % save('OrderOmega0.1.mat','SP_Traj','y','Net_Sol','tSP');
% 

hold(ax1, 'on' );
indexR = 2;
fp1=plot(ax1,tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0 0 0 0.5])
hold on 
fp2=plot(ax1,ctspanR,y(indexR,:),'--','LineWidth',3,'color','green')
xlabel(ax1,'$t \,[$s$]$','Interpreter','latex');
ylabel(ax1,'$x \,[$m$]$','Interpreter','latex');
title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$');
title(ax1,title_string,'Interpreter','latex')
hold on

grid on ;

% ax2 = axes( 'Position' ,[.2 .7 .5 .2])
hold(ax2, 'on' );
hold on 
indexR = 2;
plot(ax2,tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0 0 0 0.5])
hold on 
plot(ax2,ctspanR,y(indexR,:),'--','LineWidth',3,'color','green')
% xlabel(ax2,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax2,'$x \,[$m$]$','Interpreter','latex');


grid on ;

hold(ax3, 'on' );
hold on 
indexR = 2;
plot(ax3,tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0 0 0 0.5])
hold on 
plot(ax3,ctspanR,y(indexR,:),'--','LineWidth',3,'color','green')
% xlabel(ax3,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax3,'$x \,[$m$]$','Interpreter','latex');

% axis([0,0.8,-1.2,1.2])
grid on ;

Yrec = [Yrec,y.'];
FullTra = [FullTra,Full_Traj];

% NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));
% NMTET = [NMTET,NMTE];
    end
    store_Traj(indk+1,:) = {ctspanR,Yrec,tSP,FullTra};
end
hold(ax1, 'on' );
hold on 
fp3 =plot(ax1,alphaTOG/epsilon,XIU(:,indexR),'-','LineWidth',3,'color',[1 0 0 0.3])
hold on 
fp4 = plot(ax1,alphaTOG/epsilon,XISP(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
plot(ax1,alphaTOG/epsilon,XISM(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
legend([fp1 fp2 fp3 fp4],'Full order model','Reduced order model','Chaotic anchor trajectory','Chaotic attractors','Interpreter','latex','Location','best')

hold(ax2, 'on' );
hold on
plot(ax2,alphaTOG/epsilon,XIU(:,indexR),'-','LineWidth',3,'color',[1 0 0 0.3])
hold on 
plot(ax2,alphaTOG/epsilon,XISP(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
plot(ax2,alphaTOG/epsilon,XISM(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
% axis([0,750,0.97,1.025])

hold(ax3, 'on' );
hold on
plot(ax3,alphaTOG/epsilon,XIU(:,indexR),'-','LineWidth',3,'color',[1 0 0 0.3])
hold on 
plot(ax3,alphaTOG/epsilon,XISP(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
plot(ax3,alphaTOG/epsilon,XISM(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
% axis([0,0.8,-1.2,1.2])


% origUnits = fig.Units;
% fig.Units = fig.PaperUnits;
% fig.PaperSize = fig.Position(3:4);
% fig.Units = origUnits;
% exportgraphics(fig, 'e008xcdot_slow.pdf');

% NMTEavg = sum(NMTET)/6;

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
% print(fig, 'rail_x_slow_clarity.pdf' , '-dpdf' , '-r300' )

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'mixed_mode_error_force.pdf' , '-dpdf' , '-r300' )