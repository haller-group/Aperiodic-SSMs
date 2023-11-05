% [M,C,K,fnl,f_0,outdof,PlotFieldonDefMesh] = build_model(nElements);
% clear all
c = 0.3;
k = 1;
a = 1;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
K = [a^2*(1 + m/M)*g*4 - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];


B = eye(4);

% quadratic part
F2 = sptensor([4 4 4]);
F3 = sptensor([4 4 4 4]);
F3(3,1,1,1) = -4*g*(1+m/M);
F3(3,1,3,3) = -16*(1+m/M)*a^4;
F = {F2,F3};

DS = DynamicalSystem();
set(DS,'A',A,'B',B,'fnl',F); 
set(DS.Options,'notation','tensor')
% spectrum analysis
[V,D,W] = DS.linear_spectral_analysis();
%%
% DS set up - SSM tool

% A = V*diag(D)*inv(V);
% VI = inv(V);
% 
% V_new = [V(:,4),V(:,1),V(:,2),V(:,3)];
% D_new = [D(4),D(1),D(2),D(3)];
% VI_new = inv(V_new);
% A = V_new*diag(D_new)*VI_new;
% 
% set(DS.Options,'notation','tensor')

% w1 = sqrt(k/m);
% w2 = sqrt(3*k/m);
% ls = (-c/(2*m*w1) + 1i*sqrt(1-c^2/(4*m^2*w1^2)))*sqrt(k/m);
% lf = (-(3*c)/(2*m*w2) + 1i*sqrt(1-(3*c)^2/(4*m^2*w2^2)))*sqrt(3*k/m);

% subs2 = [1 1 1 1];
% n = 3;
% F3 = sptensor(subs2, gamma, [n,n,n,n]);
% F2 = sptensor([n,n,n]);
% fnl = {F2,F3};
% 
% 
% n = length(M);
% disp(['Number of degrees of freedom = ' num2str(n)])
% disp(['Phase space dimensionality = ' num2str(2*n)])
% 
% %%
% % DS set up - SSM tool
% hold on
% DS = DynamicalSystem();
% set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
% set(DS.Options,'Emax',6,'Nmax',10,'notation','tensor')
% [V,D,W] = DS.linear_spectral_analysis();

resonant_modes = [4 1];
Acheck = V*diag(D)*inv(V);
VI = inv(V);
% real(D(3))/real(D(1))

S = SSM(DS);
set(S.Options,'paramStyle', 'graph','reltol', 1,'notation','tensor');
set(S.Options,'notation','tensor');
                       % master mode
order = 15;                              % SSM expansion order
S.choose_E(resonant_modes);
[W0,R0] = S.compute_whisker(order);     % compute of SSM

epsilon = 0;
tf = 500;
nsteps = 8000;
q0 = 0;
q0 = [q0;conj(q0)+0.1];
z0 = reduced_to_full_traj(0,q0,W0);
traj = transient_traj_on_auto_ssm(DS, resonant_modes, W0, R0, tf, nsteps, 1:4, [], q0);


ctspan = linspace(0,500,8001);
epsilon = 0;
Force_Lorenz = @(t) 0;
GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon);
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
indexR = 1; 
plot(tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(tSP,traj.phy(:,indexR),'--','LineWidth',3,'color','red')
