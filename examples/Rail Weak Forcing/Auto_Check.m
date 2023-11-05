clear all
c = 0.3;
k = 1;
a = 0.5;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
scale = 10*a^3;%2*a^3;
K = [-2*a^2*(1 + m/M)*g*4/scale - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];


B = eye(4);

B = eye(4);
signe = 1;
% quadratic part
F2 = sptensor([4 4 4]);
F3 = sptensor([4 4 4 4]);
F3(3,1,1,1) = -4*g/scale*(1+m/M);
F3(3,1,3,3) = -16*(1+m/M)*a^4/scale^2;
F2(3,1,1) =signe*-4*g*(1+m/M)*3*a/scale;
F2(3,3,3) = signe*-16*(1+m/M)*a^5/scale^2;

F = {F2,F3};

DS = DynamicalSystem();
set(DS,'A',A,'B',B,'fnl',F); 
set(DS.Options,'notation','tensor')


% spectrum analysis
[Vs,Ds,Ws] = DS.linear_spectral_analysis();
resonant_modes = [1 2];

S = SSM(DS);
set(S.Options,'paramStyle', 'graph','reltol', 1,'notation','tensor');
set(S.Options,'notation','tensor');
                       % master mode
order = 9;                              % SSM expansion order
S.choose_E(resonant_modes);
[W0,R0] = S.compute_whisker(order);     % compute of SSM

epsilon = 0;
tf = 500;
nsteps = 8000;
q0 = 1*exp(1i*0.5);
q0 = [q0;conj(q0)];
z0 = reduced_to_full_traj(0,q0,W0);
traj = transient_traj_on_auto_ssm(DS, resonant_modes, W0, R0, tf, nsteps, 1:4, [], q0);


% ctspan = linspace(0,500,8001);
% 
% GH = @(t,y) Full_Order_Model_PM_stable(t,y,0,0,c,k,a,g,m,M,cf,signe);
% tic
% IC = traj.phy(1,:).';
% [tSP,Full_Traj] = ode45(GH,ctspan,IC);
% toc
rhosamp = [0 0.05 0.08:0.01:0.7];
plotdofs = [2 3 4]; 
figure
plot_2D_auto_SSM(0,[signe*a,signe*a*m/(M+m),0,0],W0,rhosamp,plotdofs,{'$x_c$','$\dot{x}$','$\dot{x}_c$'});
view([-70 7])
figssm = gcf;

lamdMaster = DS.spectrum.Lambda(resonant_modes);
options = struct();
options.isauto = true; 
options.isdamped = true;
options.numDigits = 4;
disp('Reduced dynamics on the 2D unstable SSM:')
y = reduced_dynamics_symbolic(lamdMaster,R0,options)


c = 0.3;
k = 1;
a = 0.3;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;


kf = 4;
scale = 10*a^3;%2*a^3;

K = [a^2*(1 + m/M)*g*4/scale - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];


B = eye(4);

% quadratic part
F2 = sptensor([4 4 4]);
F3 = sptensor([4 4 4 4]);
F3(3,1,1,1) = -4*g*(1+m/M)/scale;
F3(3,1,3,3) = -16*(1+m/M)*a^4/scale^2;
F = {F2,F3};

DS = DynamicalSystem();
set(DS,'A',A,'B',B,'fnl',F); 
set(DS.Options,'notation','tensor')
% spectrum analysis
[V,D,W] = DS.linear_spectral_analysis();
%%
% DS set up - SSM tool

A = V*diag(D)*inv(V);
VI = inv(V);

V_new = [V(:,4),V(:,1),V(:,2),V(:,3)];

D_new = [D(4),D(1),D(2),D(3)];
VI_new = inv(V_new);
A = V_new*diag(D_new)*VI_new;

tspan = linspace(0,500,100000);

l1 = D_new(1);
l2 = D_new(2);
l3 = D_new(3);
u0 = [-0.1;0.4];
y0 = compute_SSM_autonomous(V_new,VI_new,u0,a,m,M,l1,l2,l3,g,scale);

RomH = @(t,uk) rom_auto(t,uk,l1,l2,l3,V_new,VI_new,a,m,M,g,scale);
tic
IC = u0;
[tred,ROM_SOL] = ode45(RomH,tspan,IC);
toc
nssm = 20;
y_sol = compute_SSM_autonomous(V_new,VI_new,ROM_SOL.',a,m,M,l1,l2,l3,g,scale);
GH = @(t,y) Full_Order_Model_PM_vec(t,y,0,0,c,k,a,g,m,M,cf,scale);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,tspan,IC);
toc

[U1,U2] = meshgrid(linspace(-1,1,nssm));

Chart = [U1(:),U2(:)];
[U11,U21] = meshgrid(linspace(-0.05,0.05,nssm));

Chart_complex = [U11(:)+1i*U21(:),U11(:)-1i*U21(:)];
y_ssm = compute_SSM_autonomous(V_new,VI_new,Chart.',a,m,M,l1,l2,l3,g,scale);

tic
IC = [y_ssm(1,:).';y_ssm(2,:).';y_ssm(3,:).';y_ssm(4,:).'];
[tSP,Full_SSM] = ode45(GH,tspan,IC);
toc



 Z1 = reshape(y_ssm(2,:),nssm,nssm);
     Z2 = reshape(y_ssm(3,:),nssm,nssm);
     Z3 = reshape(y_ssm(4,:),nssm,nssm);
planeUS = real(V_new*[Chart.';Chart.'*0]);
planeS = real(Vs*[Chart_complex.';Chart_complex.'*0]);

figure 
h = surf(Z1,Z2,Z3)
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     hold on
     plot3(Full_Traj(:,2),Full_Traj(:,3),Full_Traj(:,4),'-','LineWidth',3,'color','blue')
     hold on 
     plot3(y_sol(2,:),y_sol(3,:),y_sol(4,:),'--','LineWidth',3,'color','red')
     xlabel('$x \,[$m$]$','Interpreter','latex');
ylabel('$x_c \,[$m$]$','Interpreter','latex');
zlabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');

figure
plot3(Full_Traj(:,1),Full_Traj(:,2),Full_Traj(:,3),'-','LineWidth',3,'color','blue')
xlabel('$x \,[$m$]$','Interpreter','latex');
ylabel('$x_c \,[$m$]$','Interpreter','latex');
zlabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');


figure 
 Z1 = reshape(y_ssm(4,:),nssm,nssm);
  h = surf(U1,U2,Z1)
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     hold on

figure 

Z1 = reshape(y_ssm(1,:),nssm,nssm);
     Z2 = reshape(y_ssm(3,:),nssm,nssm);
     Z3 = reshape(y_ssm(4,:),nssm,nssm);
h = surf(Z1,Z2,Z3)
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;


figure
indexR=1;
plot(tspan,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
plot(tspan,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$x \,[$m$]$','Interpreter','latex');

figure
indexR=4;
plot(tspan,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
plot(tspan,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$x_c \,[$m$]$','Interpreter','latex');



figure 
plot(ROM_SOL(:,1),ROM_SOL(:,2))


%% video 
 Np = max(size(Full_SSM(1,:)))/4;
Data_matrix1 = Full_SSM(:,1:Np);
Data_matrix2 = Data_matrix1(:);
Data_matrix = [];
Data_matrix = [Data_matrix;Data_matrix2.'];

Data_matrix1 = Full_SSM(:,Np+1:2*Np);
Data_matrix2 = Data_matrix1(:);

Data_matrix = [Data_matrix;Data_matrix2.'];

Data_matrix1 = Full_SSM(:,2*Np+1:3*Np);
Data_matrix2 = Data_matrix1(:);

Data_matrix = [Data_matrix;Data_matrix2.'];

Data_matrix1 = Full_SSM(:,3*Np+1:4*Np);
Data_matrix2 = Data_matrix1(:);

Data_matrix = [Data_matrix;Data_matrix2.'];

[u,s,v] = svds(Data_matrix, 2);
V = (s\u'./max(abs(v'),[],2))'; % tangent space
iMmap = @(y) V'*y;
Xi = iMmap(Data_matrix);

% figure
% for ind = 1:50
% plot(Xi(1,1+max(size(Xi(1,:)))/400*(ind-1):ind*max(size(Xi(1,:))/400)),Xi(2,1+max(size(Xi(1,:)))/400*(ind-1):ind*max(size(Xi(1,:))/400)),'LineWidth',3)
% hold on 
% end

POD_Plane = V*Chart.';

figure 
redu = VI_new*Data_matrix;
plot(real(redu(1,1:max(size(redu(1,:))/400))),real(redu(2,1:max(size(Xi(2,:))/400))),'LineWidth',3)

figure 
redu = inv(Vs)*Data_matrix;
plot(real(redu(1,1:max(size(redu(1,:))/400))),imag(redu(1,1:max(size(Xi(2,:))/400))),'LineWidth',3)

reshape(A', [], 1)
hFig = figure('DefaultAxesFontSize',20);                       % Bring up new figure
imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
hFig.WindowState = 'maximized';
pause(0.5)
delta = 0.5;
tl =  50;
clear movieVector
ind = 1;
for j = 1:1%100%(max(size(t_sol))-tl)
i = j+0;
     clf
% subplot(2, 2, [1 3])
    Np = max(size(Full_SSM(1,:)))/4;
    ZZ = [Full_SSM(j,1:Np);Full_SSM(j,Np+1:2*Np); Full_SSM(j,2*Np+1:3*Np); Full_SSM(j,3*Np+1:4*Np)];
 plotdof = [2,3,4];
     Z1 = reshape(ZZ(plotdof(1),:),sqrt(Np),sqrt(Np));
     Z2 = reshape(ZZ(plotdof(2),:),sqrt(Np),sqrt(Np));
     Z3 = reshape(ZZ(plotdof(3),:),sqrt(Np),sqrt(Np));

     h = surf(Z1,Z2,Z3)
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     hold on
     grid on 

     planeUS1 = reshape(planeUS(plotdof(1),:),sqrt(Np),sqrt(Np));
     planeUS2 = reshape(planeUS(plotdof(2),:),sqrt(Np),sqrt(Np));
     planeUS3 = reshape(planeUS(plotdof(3),:),sqrt(Np),sqrt(Np));

     h = surf(planeUS1,planeUS2,planeUS3)
     h.EdgeColor = 'none';
     h.FaceColor = 'red';
     h.FaceAlpha = 0.3;
     hold on
     grid on 

      hold on

      podUS1 = reshape(POD_Plane(plotdof(1),:),sqrt(Np),sqrt(Np));
     podUS2 = reshape(POD_Plane(plotdof(2),:),sqrt(Np),sqrt(Np));
     podUS3 = reshape(POD_Plane(plotdof(3),:),sqrt(Np),sqrt(Np));

     h = surf(podUS1,podUS2,podUS3)
     h.EdgeColor = 'none';
     h.FaceColor = 'green';
     h.FaceAlpha = 0.3;
     hold on
     grid on 

      hold on
     plot3(Full_Traj(:,plotdof(1)),Full_Traj(:,plotdof(2)),Full_Traj(:,plotdof(3)),'-','LineWidth',3,'color','blue')
      planeS1 = reshape(planeS(plotdof(1),:),sqrt(Np),sqrt(Np));
     planeS2 = reshape(planeS(plotdof(2),:),sqrt(Np),sqrt(Np));
     planeS3 = reshape(planeS(plotdof(3),:),sqrt(Np),sqrt(Np));

     h = surf(planeS1,planeS2,planeS3)
     h.EdgeColor = 'none';
     h.FaceColor = 'blue';
     h.FaceAlpha = 0.3;
     hold on
     grid on 

     hold on 
     rhosamp = linspace(0,0.02,20);
plotdofs = [2 3 4]; 
plot_2D_auto_SSM(0,[signe*a,signe*a*m/(M+m),0,0],W0,rhosamp,plotdofs,{'$x_c$','$\dot{x}$','$\dot{x}_c$'});
view([-70 7])
figssm = gcf;

%      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
     starto = 1;
     
%      plot3(SolutionZ(plotdofs(1),starto),SolutionZ(plotdofs(2),starto),SolutionZ(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%     hold on
%     plot3(SP_TrajZ(plotdofs(1),j:j+tl),SP_TrajZ(plotdofs(2),j:j+tl),SP_TrajZ(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color','black');
%      hold on
%      plot3(SP_TrajZ(plotdofs(1),j:j+tl),SP_TrajZ(plotdofs(2),j:j+tl),SP_TrajZ(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color','black');
%      hold on
%      plot3(shift_trajZ(plotdofs(1),j:j+tl),shift_trajZ(plotdofs(2),j:j+tl),shift_trajZ(plotdofs(3),j:j+tl),'--','LineWidth',3,'Color','red');
%      hold on
%      plot3(SolutionZ(plotdofs(1),j:j+tl),SolutionZ(plotdofs(2),j:j+tl),SolutionZ(plotdofs(3),j:j+tl),'-','LineWidth',3,'color',[0,1,0,0.3])
%      hold on
%      plot3(shift_trajZ(plotdofs(1),j+tl),shift_trajZ(plotdofs(2),j+tl),shift_trajZ(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','red');
%      hold on
%      plot3(SP_TrajZ(plotdofs(1),j+tl),SP_TrajZ(plotdofs(2),j+tl),SP_TrajZ(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','black');
%      hold on
%       
%     xlabel('$Re(u_1+\xi_{u_1})$','Interpreter','latex');
%      ylabel('$Imag(u_1 + \xi_{u_1})$','Interpreter','latex');
%      zlabel('$Re(v_2 + \xi_{v_2})$','Interpreter','latex');
% %      title('Modal Coordinates SSM for \epsilon = 0.1 and |\epsilon F_{max}| = 10')
%     legend('SSM - W_{\epsilon}(x_{\epsilon}(t)) O(\epsilon^{k_1} u^{k_2}, k_1 + k_2 = 3)','FOM','ROM','O(\epsilon^3) Anchor Trajectory')
% %     axis([-0.5,0.5,-0.5,0.5,-0.02,0.02])
% %      view(-116,17)
%     grid on
%        % daspect([1,1,1])
%  axis([-0.2,0.2,-0.2,0.2,-0.01,0.01])
% set(gcf,'color','white')
%     figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

end