clear all
c = 0.3;
k = 1;
a = 0.3;
g = 9.8;
m = 1;
M = 4;
cf = 0.03;
K = [-2*a^2*(1 + m/M)*g*4 - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];


B = eye(4);
signe = -1;
% quadratic part
F2 = sptensor([4 4 4]);
F3 = sptensor([4 4 4 4]);
F3(3,1,1,1) = -4*g*(1+m/M);
F3(3,1,3,3) = -16*(1+m/M)*a^4;
F2(3,1,1) =signe*-4*g*(1+m/M)*3*a;
F2(3,3,3) = signe*-16*(1+m/M)*a^5;

F = {F2,F3};



DS = DynamicalSystem();
set(DS,'A',A,'B',B,'fnl',F); 
set(DS.Options,'notation','tensor')


% spectrum analysis
[Vs,D,W] = DS.linear_spectral_analysis();
%%
% DS set up - SSM tool

A = V*diag(D)*inv(V);
VI = inv(V);



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
q0 = 100*exp(1i*0.5);
q0 = [q0;conj(q0)];
z0 = reduced_to_full_traj(0,q0,W0);
traj = transient_traj_on_auto_ssm(DS, resonant_modes, W0, R0, tf, nsteps, 1:4, [], q0);


ctspan = linspace(0,500,8001);

GH = @(t,y) Full_Order_Model_PM_stable(t,y,0,0,c,k,a,g,m,M,cf,signe);
tic
IC = traj.phy(1,:).';
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc


lamdMaster = DS.spectrum.Lambda(resonant_modes);
options = struct();
options.isauto = true; 
options.isdamped = true;
options.numDigits = 4;
disp('Reduced dynamics on the 2D unstable SSM:')
y = reduced_dynamics_symbolic(lamdMaster,R0,options)

figure 
indexR = 4; 
plot(tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(tSP,traj.phy(:,indexR),'--','LineWidth',3,'color','red')

figure 
plot3(Full_Traj(:,1),Full_Traj(:,2),Full_Traj(:,3),'-','LineWidth',3,'color','black');
hold on
plot3(traj.phy(:,1),traj.phy(:,2),traj.phy(:,3),'-','LineWidth',3,'color','red');

%%

c = 0.3;
k = 1;
a = 0.3;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
kf = 4;
scale = 1;

ctspan = linspace(0,500,12001);

GH = @(t,y) Full_Order_Model_PM(t,y,0,0,c,k,a,g,m,M,cf,kf,scale);
tic
IC = [0;0;0;0.1];
[tSP,Full_TrajU] = ode45(GH,ctspan,IC);
toc
fig = figure 
plot3(Full_TrajU(:,1),Full_TrajU(:,2),Full_TrajU(:,4),'-','LineWidth',3,'color','black');
tic
IC = [0;0;0;-0.1];
[tSP,Full_TrajU] = ode45(GH,ctspan,IC);
toc
hold on 
plot3(Full_TrajU(:,1),Full_TrajU(:,2),Full_TrajU(:,4),'-','LineWidth',3,'color',[0.8500 0.3250 0.0980]);

hold on 
tic
IC =  0.1*real(V(:,1));
[tSP,Full_TrajU] = ode45(GH,ctspan,IC);
toc
hold on 
plot3(Full_TrajU(:,1),Full_TrajU(:,2),Full_TrajU(:,4),'-','LineWidth',3,'color','green');

hold on 
tic
IC = -0.01*real(V(:,2));
[tSP,Full_TrajU] = ode45(GH,ctspan,IC);
toc
hold on 
plot3(Full_TrajU(:,1),Full_TrajU(:,2),Full_TrajU(:,4),'-','LineWidth',3,'color','green');

hold on 
plot3(0,0,0,'.','MarkerSize',30,'color','red')
hold on 
plot3(a,a*m/(m+M),0,'.','MarkerSize',30,'color','blue')
hold on 
plot3(-a,-a*m/(m+M),0,'.','MarkerSize',30,'color','blue')
xlabel('$x$','Interpreter','latex')
zlabel('$\dot{x}_c$','Interpreter','latex')
%%
% axis([-5 5 -5 5 -5 5])
figure
indexR = 1;
plot(tSP,Full_TrajU(:,indexR),'-','LineWidth',3,'color','black')




