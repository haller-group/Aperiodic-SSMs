
% [M,C,K,fnl,f_0,outdof,PlotFieldonDefMesh] = build_model(nElements);
hold on 
m = 1;
k = 1;
c = 0.3;
gamma = 0.5;
M = [m,0;0,m];
C = [c,-c;-c,2*c];
K = [2*k,-k;-k,2*k];
n = length(M);

subs2 = [1 1 1 1];

F3 = sptensor(subs2, gamma, [n,n,n,n]);
F2 = sptensor([n,n,n]);
fnl = {F2,F3};


n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])

%%
% DS set up - SSM tool
hold on
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
real(D(3))/real(D(1))
S = SSM(DS);
resonant_modes = [1 2];
set(S.Options, 'reltol', 1,'notation','multiindex');
order = 7;
S.choose_E(resonant_modes)
[W0,R0] = S.compute_whisker(order);

lamdMaster = DS.spectrum.Lambda(resonant_modes);
options = struct();
options.isauto = true; 
options.isdamped = true;
options.numDigits = 6;
disp('Reduced dynamics on the 2D stable SSM:')
y = reduced_dynamics_symbolic(lamdMaster,R0,options)


rhosamp = linspace(0,1,201);
     plotdofs = [1 2 4]; 
      
%      W0 = WTF(target(Omega*traj.time(i+tl)));
  
coord_change_vec2 = 0;
coord_change_vec = [0;0;0];
[zdof1,zdof2,zdof3]=plot_2D_auto_SSM(coord_change_vec2,coord_change_vec,W0,rhosamp,plotdofs,{'$q_1$ m','$q_2$ m','$\dot{q}_2$ m/s'});
   

%%
% Define forcing - transverse forcing along the beam length. 
% Term this mildly periodic - for long times it will look periodic (backward and forward)
% F_{alpha} = [atan(alpha)/(pi/2) +Asin(alpha)] ,max_force = 1+A, alpha = Omega*time + t_0 (offset) 
% Twist plan - F_{alpha} = A*atan(alpha)/(pi/2) , F_{t} = epsilon*sin(w*t),
% F_{net} = F_{alpha} + F_{t}
% For now tuning parameters based on plots
% Expectation - transistion from one limit cycle to the other (orbit)
Omega = 0.1;
% f_tspan = linspace(0,16*pi/(Omega),1000);
% A = 0.5;
% t0 = -25;
% load_vector = zeros(n,1);
% transverse_coor = 2:3:(n-1);
% load_vector(transverse_coor,1) = 1;
% F_alpha = epsilon*(2*atan(Omega*f_tspan+t0)./(pi/2) + A*sin(Omega*f_tspan+t0));
% Load_full = load_vector.*F_alpha;
% figure
% plot(f_tspan,Load_full(2,:))


%Lorenz force 
epsilon = 1;
tspan = linspace(0,60,1000);
loren = @(t,y) lorenz_3D(t,y);
IC = [0.8;0.3;0.2];
[tk,yy] = ode45(loren,tspan,IC);
F_a = yy(20:end,1)/(max(yy(20:end,1)));
F_a_new=2*epsilon*F_a;
tspan = tspan(20:end)-tspan(20);
Force_Lorenz = @(xq) interp1(tspan(:),F_a_new(:),xq,'spline');
sigma = 10;
rho = 28;
beta = 8/3;
% DF_a_new = 2*epsilon.*(sigma*(yy(:,1).*(rho - yy(:,3))-yy(:,2)-sigma.*(yy(:,2)-yy(:,1)))/(max(yy(:,1))));

[dXdFPa, XtruncPa,F_truncPa] = ftd(F_a_new.', tspan);
 
DForce_Lorenz = @(xq) interp1(F_truncPa(:),dXdFPa(:),xq,'spline');



% Load_full_a = zeros(n,1000);
% Load_full_a(2:3:n,:) = repmat(F_a_new.',nElements,1);
% 
% figure
% plot(tk,F_a_new,'.')


%% Quasi Static Response for transverse loads
tic
[F, lambda, V, G, DG] = functionFromTensors(M, C, K, fnl);
toc
disp('Starting Quasi-Static response estimator') % 33 minutes roughly to compute the estimator
% % tic
% % IC = getStaticResponse(K, M, F, Load_full_a, 1, PlotFieldonDefMesh);
% % toc
% % ICL = IC(1:n,:);
% % save('Quasi_Static_Lorenz.mat','ICL')
% % load('Quasi_Static.mat','IC')
% 

 N_grid = 1000;
% Grid_F = linspace(-max(F_a_new),max(F_a_new),N_grid);
% WT = cell(1,N_grid);
% RT = cell(1,N_grid);
% DJ = cell(1,N_grid);
% ICTT = cell(1,N_grid);
% spec_ratio = [];
% backbonecurve = cell(1,N_grid);
% for i = 1:N_grid
% load_vector = zeros(n,1);
% load_vector(2:3:n,1) = Grid_F(i);
% IC = getStaticResponse(K, M, F, load_vector, 0, 0);
% ICTT{1,i} = IC(1:n,1);
% [K_shift,fnl_shift] = build_model_shift(K,fnl,IC(1:n,1));
% DS = DynamicalSystem();
% set(DS,'M',M,'C',C,'K',K_shift,'fnl',fnl_shift);
% set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% [V,D,W] = DS.linear_spectral_analysis();
% V = reorient_vec(V);
% W = reorient_vec(W);
% DS.spectrum.V = V;
% DS.spectrum.W = W;
% 
% spec_ratio = [spec_ratio,real(D(3))/real(D(1))];
% S = SSM(DS);
% set(S.Options, 'reltol', 1,'notation','multiindex');
   resonant_modes = [1 2]; % choose master spectral subspace
% order = 5;                  % SSM expansion order
% S.choose_E(resonant_modes)
% [W0,R0] = S.compute_whisker(order);
% WT{1,i} = W0;
% RT{1,i} = R0;
% lamdMaster = DS.spectrum.Lambda(resonant_modes);
% options = struct();
% options.isauto = true; 
% options.isdamped = true;
% options.numDigits = 6;
% disp('Reduced dynamics on the 2D stable SSM:')
% y = reduced_dynamics_symbolic(lamdMaster,R0,options)
% DJ{1,i} = lamdMaster;
% 
% % hj = matlabFunction(y(2));
% % backbonecurve{1,i} = hj;
% end    
% 
%  save('RTSP.mat','RT','R0','Grid_F');
%  save('WTSP.mat','WT','W0','Grid_F');
%  save('DSP.mat','DJ','Grid_F');
%  save('ICSP.mat','ICTT','Grid_F');


%  lambda = [];
%  for i = 1:1000
%     lambda = [lambda,real(DJ{1,i}(1,1))];
%  end    


% X_s = IC(1:n,1);
% [K_shift,fnl_shift] = build_model_shift(K,fnl,X_s);
% DS = DynamicalSystem();
% set(DS,'M',M,'C',C,'K',K_shift,'fnl',fnl_shift);
% set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% [V,D1,W] = DS.linear_spectral_analysis();

%%

% load('RTbeam4.mat','RT','R0','Grid_F');
% load('WTbeam4.mat','WT','W0','Grid_F');
% load('Dbeam4.mat','DJ','Grid_F');
% load('ICbeam4.mat','ICTT','Grid_F');

load('RTSP.mat','RT','R0','Grid_F');
load('WTSP.mat','WT','W0','Grid_F');
load('DSP.mat','DJ','Grid_F');
load('ICSP.mat','ICTT','Grid_F');

RF = grid_SSM_interp1(RT,Grid_F,N_grid);
WF = grid_SSM_interp1(WT,Grid_F,N_grid);
IF = grid_SSM_interp1_nc(ICTT,Grid_F,N_grid);

RTF = @(xa) convert_fit_to_value_x(RF,xa,R0);
WTF = @(xa) convert_fit_to_value_x(WF,xa,W0);

ICTTF = @(xa) convert_fit_to_value_x_na(IF,xa,zeros(n,1));

Load_Check = cell2mat(ICTT);
[dXdF, Xtrunc,F_trunc] = ftd(Load_Check, Grid_F);
Boad_Check = cell(1,max(size(F_trunc)));
for ind = 1:max(size(F_trunc))
Boad_Check{1,ind} = dXdF(:,ind);
end    
DIF = grid_SSM_interp1_nc(Boad_Check,F_trunc,max(size(F_trunc)));
DICTTF = @(xa) convert_fit_to_value_x_na(DIF,xa,zeros(n,1));
% % 
Omega = 0.01;
q0 = 0.6*exp(1i*0.5);
q0 = [q0;conj(q0)];
target = @(alpha) {Force_Lorenz(alpha)};
z0 = reduced_to_full_traj_adiabatic_sgrid(0,q0,WTF,Omega,target);
tf = 58.8589/Omega;
nsteps = 20000;
outdof = [29];

tic
traj3 = transient_traj_on_auto_ssm_adiabatic_sgrid(DS, resonant_modes, WTF, RTF, Omega,target,tf, nsteps, outdof, [], q0);
toc
% % % 
 Gh = @(t,x) F(t,x)+[0;0;0;Force_Lorenz(Omega*t)];
% % % % % 
tic
ctspan = linspace(0,tf,nsteps + 1);
IC = (z0 +[ICTTF(target(0));DICTTF(target(0))*Omega*DForce_Lorenz(0)]);
[tSP,SP_Traj] = ode45(Gh,ctspan,IC);
toc


tic
ctspan = linspace(0,tf,20000);
IC = (z0 +[ICTTF(target(0));DICTTF(target(0))*Omega*DForce_Lorenz(0)])+[4;4;4;4];
[tSP_Off,SP_Traj_Off] = ode45(Gh,ctspan,IC);
toc


 shift_traj = shift_to_target(DICTTF,ICTTF,traj3,Omega,target,DForce_Lorenz);
% % % 
traj_actual = traj3.phy;
for ind = 1:max(size(traj3.time))
   traj_actual(ind,:) =[ICTTF(target(Omega*traj3.time(ind)));DICTTF(target(Omega*traj3.time(ind)))*Omega*DForce_Lorenz(Omega*traj3.time(ind))].';
end
% % 

NMTE = sum(sqrt(sum((shift_traj - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(tSP)))
RelErrorO4 = (sqrt(sum((shift_traj - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));

figure
plot(alpha,RelErrorO2,'-','LineWidth',3)
hold on
plot(alpha,RelErrorO3,'-','LineWidth',3)
hold on
plot(alpha,RelErrorO1,'-','LineWidth',3)
hold on 
plot(alpha,RelErrorO4,'-','LineWidth',3)

figure 

plot(tSP,SP_Traj(:,1),'-','LineWidth',3,'color','black');
hold on
plot(traj3.time,shift_traj(:,1),'--','LineWidth',3,'color','red');
hold on 
plot(traj3.time,traj_actual(:,1),'-','LineWidth',3,'color',[0 1 0 0.3]);
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$q_1 \,[$m$]$','Interpreter','latex');

figure 
subplot(2,2,1)
plot(tSP,SP_Traj(:,1),'-','LineWidth',3,'color','black');
hold on
plot(traj3.time,shift_traj(:,1),'--','LineWidth',3,'color','red');
hold on 
plot(traj3.time,traj_actual(:,1),'-','LineWidth',3,'color',[0 1 0 0.3]);
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$q_1 \,[$m$]$','Interpreter','latex');

legend('FOM','ROM','Target')
 

subplot(2,2,2)
plot(tSP,SP_Traj(:,2),'-','LineWidth',3,'color','black');
hold on
plot(traj3.time,shift_traj(:,2),'--','LineWidth',3,'color','red');
hold on 
plot(traj3.time,traj_actual(:,2),'-','LineWidth',3,'color',[0 1 0 0.3]);
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$q_2 \,[$m$]$','Interpreter','latex');




subplot(2,2,3)

plot(tSP,SP_Traj(:,3),'-','LineWidth',3,'color','black');
hold on
plot(traj3.time,shift_traj(:,3),'--','LineWidth',3,'color','red');
hold on 
plot(traj3.time,traj_actual(:,3),'-','LineWidth',3,'color',[0 1 0 0.3]);
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$\dot{q}_1 \,[$m/s$]$','Interpreter','latex');




subplot(2,2,4)
plot(tSP,SP_Traj(:,4),'-','LineWidth',3,'color','black');
hold on
plot(traj3.time,shift_traj(:,4),'--','LineWidth',3,'color','red');
hold on 
plot(traj3.time,traj_actual(:,4),'-','LineWidth',3,'color',[0 1 0 0.3]);


xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$\dot{q}_2 \,[$m/s$]$','Interpreter','latex');



sgtitle('\Omega = 0.01')

% % % 
% % 
% % figure 
% % plot(linspace(0,60/(0.1),3001),xi0,'-*','LineWidth',3,'color','black');
% % hold on
% % plot(linspace(0,60/(0.05),3001),xi1,'-o','LineWidth',3,'color','red');
% % hold on 
% % plot(linspace(0,60/(0.01),3001),xi2,'-square','color','blue');
% % hold on
% % plot(linspace(0,60/(0.001),3001),xi3,'-+','color','green');
% % 
% % xlabel('$t \,[$s$]$','Interpreter','latex');
% % ylabel('$|\xi|$','Interpreter','latex');
% % 
% % axis([0 600 0 0.9])
% % 
% % 
% % figure 
% % plot(linspace(0,60,3001),xi0,'-*','LineWidth',3,'color','black');
% % hold on
% % plot(linspace(0,60,3001),xi1,'-o','LineWidth',3,'color','red');
% % hold on 
% % plot(linspace(0,60,3001),xi2,'-square','color','blue');
% % hold on
% % plot(linspace(0,60,3001),xi3,'-+','color','green');
% % 
% % xlabel('$\alpha$','Interpreter','latex');
% % ylabel('$|\xi|$','Interpreter','latex');
% % 
% % axis([0 60 0 0.9])

% % % 
% % % tic
% % % tspan = linspace(0,10,1000);
% % % IC = (z0*0 +[ICTTF(target(0));0*DICTTF(target(0))*Omega*DForce_Lorenz(0)]);
% % % [tBeam,Beam_SSM_Traj] = ode15s(Gh,tspan,IC);
% % % toc
% % % 
% % % tic
% % % tspan = linspace(0,tf,nsteps);
% % % Off_epsilon = 0.02;
% % % IC = IC + Off_epsilon;
% % % [Off_tBeam,Off_Beam_SSM_Traj] = ode45(Gh,tspan,IC);
% % % toc
% % % 







