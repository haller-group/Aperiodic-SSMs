

fig = figure
fontsize(fig, 20, "points")
set(gcf,'color','w');
box on
grid on ;
axis([0,150,-0.5,1])
ax1 = gca;
ax2 = axes( 'Position' ,[.2 .7 .35 .2]); 
box on
grid on ;
axis([0,400,0.278,0.322])
ax3 = axes( 'Position' ,[.65 .7 .2 .2]); 
box on
grid on ;
axis([60,62,-0.45,0.45])

NMTET = [];
% z0 = [0.1*exp(1i*0),0.1*exp(1i*pi/6),0.1*exp(1i*pi/3),0.1*exp(1i*pi/2),0.1*exp(2*1i*pi/3),0.1*exp(5*1i*pi/6),0.1*exp(6*1i*pi/6),0.1*exp(7*1i*pi/6),0.1*exp(8*1i*pi/6),0.1*exp(9*1i*pi/6),0.1*exp(10*1i*pi/6),0.1*exp(11*1i*pi/6)];

b = 0.3;
z0 = [ b+1i*0, ...
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

epsilon = 1;
ctspan = linspace(0,500,80000);
Unstable_Solution = epsilon*xi_full_eval(ctspan) + epsilon^3*xi_full_eval3(ctspan);

Stable_solution1 = epsilon*XI_hyperbolic_solution_a1{1,1}(ctspan) + epsilon^2*XI_hyperbolic_solution_a1{1,2}(ctspan) + epsilon^3*XI_hyperbolic_solution_a1{1,3}(ctspan);
Stable_solution2 = epsilon*XI_hyperbolic_solution_a2{1,1}(ctspan) + epsilon^2*XI_hyperbolic_solution_a2{1,2}(ctspan) + epsilon^3*XI_hyperbolic_solution_a2{1,3}(ctspan);
Stable_solution1(:,1) = Stable_solution1(:,1)+a;
Stable_solution1(:,2) = Stable_solution1(:,2)+a*m/(M+m);

Stable_solution2(:,1) = Stable_solution2(:,1)-a;
Stable_solution2(:,2) = Stable_solution2(:,2)-a*m/(M+m);




to1=find(ctspan>=20);
to2=find(ctspan>=40);
to3=find(ctspan>=60);

to = [1,to1(1),to2(1),to3(1)];



store_Traj = cell(5,4);
store_Traj(1,:) = {ctspan,Unstable_Solution.',Stable_solution1.',Stable_solution2.'};
for indk = 1:4
Yrec = [];
FullTra = [];
    for ind = 1:16



stop_index = find(ctspan>=1.5);

ctspanR = ctspan(to(indk):to(indk)+stop_index(1));
ctspanF = ctspan(to(indk):end);

u0 = [real(z0(ind));imag(z0(ind))];

[y0,modal0] = compute_SSM_phy(XI_hyperbolic_solution,ctspanR(1),V_new,VI_new,u0,h030,h300,h120,h210,h111,h201,h021,h012,h102,epsilon);

RomH = @(t,uk) rom_temp_model_noscale(t,uk,l1,l2,V_new,VI_new,epsilon,xi_sol1,xi_sol3,h300,h030,h120,h210,Sc,a,g,m,M);
tic
IC = u0;
[tred,ROM_SOL] = ode45(RomH,ctspanR,IC);
toc

[y_sol,modal_sol] = compute_SSM_phy(XI_hyperbolic_solution,tred.',V_new,VI_new,ROM_SOL.',h030,h300,h120,h210,h111,h201,h021,h012,h102,epsilon);

GH = @(t,y) Full_Order_Model_PM_vec(t,y,Force_Lorenz,epsilon,c,k,a,g,m,M,cf,Sc);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspanF,IC);
toc


% NMTE = sum(sqrt(sum((y_sol.' - Full_Traj).^2,2)))/(sqrt(max(sum((Full_Traj).^2,2)))*max(size(ctspan)))

hold(ax1, 'on' );
indexR=1;
fp1 = plot(ax1,ctspanF,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
fp2 = plot(ax1,ctspanR,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
xlabel(ax1,'$t \,[$s$]$','Interpreter','latex');
ylabel(ax1,'$x \,[$m$]$','Interpreter','latex');
title_string = strcat('max$(|F_{ext}|) = ',num2str(epsilon*(M+m)),'$ [$\mathbf{N}$], expansion order $N = 3$');
title(ax1,title_string,'Interpreter','latex')

hold on

grid on ;


hold(ax2, 'on' );
hold on 
plot(ax2,ctspanF,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
plot(ax2,ctspanR,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
% xlabel(ax2,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax2,'$x \,[$m$]$','Interpreter','latex');


grid on ;

hold(ax3, 'on' );
hold on 
plot(ax3,ctspanF,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
plot(ax3,ctspanR,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
% xlabel(ax3,'$t \,[$s$]$','Interpreter','latex');
% ylabel(ax3,'$x \,[$m$]$','Interpreter','latex');

% axis([0,0.8,-1.2,1.2])
grid on ;


% Yrec = [Yrec;modal_sol];
% FullTra = [FullTra;VI_new*Full_Traj.'];

Yrec = [Yrec;y_sol];
FullTra = [FullTra;Full_Traj.'];


    end
    store_Traj(indk+1,:) = {ctspanR,Yrec,tSP,FullTra};
end
hold(ax1, 'on' );
hold on 
fp3 = plot(ax1,ctspan,Unstable_Solution(:,indexR),'-','LineWidth',3,'color',[1,0,0,0.5]);
hold on 
fp4 = plot(ax1,ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on 
plot(ax1,ctspan,Stable_solution2(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on
legend([fp1 fp2 fp3 fp4],'Full order model','Reduced order model','Chaotic anchor trajectory','Chaotic attractors','Interpreter','latex')


hold(ax2, 'on' );
hold on 
plot(ax2,ctspan,Unstable_Solution(:,indexR),'-','LineWidth',3,'color',[1,0,0,0.5]);
hold on 
plot(ax2,ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on 
plot(ax2,ctspan,Stable_solution2(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on

% axis([0,750,0.97,1.025])

hold(ax3, 'on' );
hold on 
plot(ax3,ctspan,Unstable_Solution(:,indexR),'-','LineWidth',3,'color',[1,0,0,0.5]);
hold on 
plot(ax3,ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on 
plot(ax3,ctspan,Stable_solution2(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.5]);
hold on

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
% print(fig, 'rail_x_weak.pdf' , '-dpdf' , '-r300' )
% 
% print(fig, 'rail_force_ratio.pdf' , '-dpdf' , '-r300' )

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'forcing_weak.pdf' , '-dpdf' , '-r300' )