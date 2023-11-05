% clear all
clear all
c = 0.3;
k = 1;
a = 0.5;
g = 9.8;
m = 1;
Mf = 4;
cf = 0.3;
scale = 10*a^3;%2*a^3;
Sc = scale;
K = -[-k,k*m/(m+Mf);k*m/(m+Mf),-8*a^2*m*g/scale-k*m^2/(m+Mf)^2];
C = -[-c,c*m/(m+Mf);c*m/(m+Mf),-cf-c*m^2/(m+Mf)^2];
M = [m+Mf,0;0,m*Mf/(m+Mf)];


sign = 1;
n = 2;
F3 = sptensor([n,n,n,n]);
F2 = sptensor([n,n,n]);
F3(2,2,2,2) = 4*g*m/scale;
F2(2,2,2) = sign*4*g*m*a*3/scale;
% F3(2,2,4,4) = 16*m*a^4;

fnl = {F2,F3};


n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])


%%
% DS set up - SSM tool
hold on
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',6,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
D
S = SSM(DS);
resonant_modes = [1 2];
Acheck = V*diag(D)*inv(V);
VI = inv(V);

tic
[F, lambda, V, G, DG] = functionFromTensors(M, C, K, fnl);
toc


%Lorenz force 
%%

L = -20;
N_Max = 100000;
tspan = linspace(L,20,N_Max);
loren = @(t,y) lorenz_3D(t,y);
IC = [0.8;0.3;0.2];
[tk,yy] = ode45(loren,tspan,IC);

sigma = 10;
rho = 28;
beta = 8/3;

% F_a = -(sigma*rho*yy(:,1) - sigma*yy(:,1).*yy(:,3) - sigma*yy(:,2) - sigma^2*yy(:,2) + sigma^2*yy(:,1));
scaling = 3;
F_a_new = scaling*yy(:,1)/max(abs(yy(:,1)));


Force_Lorenz = griddedInterpolant(tspan(:),F_a_new(:),'spline');



figure
plot(tk,F_a_new,'Color','blue')
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$F(t) \,[$N$]$','Interpreter','latex');
% axis([0 300 -1 1])

pos_t = find(tk>=0);
figure
plot(tk(pos_t),F_a_new(pos_t),'Color','red')
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$F(t) \,[$N$]$','Interpreter','latex');
axis([0 10 -10 10])

[Vcheck,Dcheck] = eig(full(Acheck));
dsa = diag(Dcheck);
[d,ind]=sort(real(dsa));
Ds=Dcheck(ind,ind);         %sorted eigenvalue matrix
Vs=Vcheck(:,ind); 
V0 = Vs;

 N_grid = 1000;
Grid_F = linspace(-max(F_a_new),max(F_a_new),N_grid);
VT = cell(1,N_grid);
DT = cell(1,N_grid);
ICTT = cell(1,N_grid);
AT = cell(1,N_grid);
for i = 1:N_grid
load_vector = [Grid_F(i);0];
IC = getStaticResponse(K, M, F, load_vector, 0, 0);
ICTT{1,i} = IC(1:n,1);
[K_shift,fnl_shift] = build_model_shift(K,fnl,IC(1:n,1));
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K_shift,'fnl',fnl_shift);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
AT{1,i} = DS.BinvA;

DT{1,i} = D;

% VC = reorient_complex(V,V0);

[Vcheck,Dcheck] = eig(full(AT{1,i}));
dsa = diag(Dcheck);
[d,ind]=sort(real(dsa));
Ds=Dcheck(ind,ind);         %sorted eigenvalue matrix
Vs=Vcheck(:,ind); 

 VC = reorient_complex(Vs,V0);
% V0 = VC;
% [V,DA] = eigs(AT{1,i});
% dsa = diag(DA);
% [d,ind]=sort(real(-dsa));
% Ds=DA(ind,ind)          %sorted eigenvalue matrix
% Vs=V(:,ind); 
VT{1,i} = VC;

% % spec_ratio = [spec_ratio,real(D(3))/real(D(1))];
% % S = SSM(DS);
% % set(S.Options, 'reltol', 1,'notation','multiindex');
% %    resonant_modes = [1 2]; % choose master spectral subspace
% % order = 5;                  % SSM expansion order
% % S.choose_E(resonant_modes)
% % [W0,R0] = S.compute_whisker(order);
% % WT{1,i} = W0;
% % RT{1,i} = R0;
% % lamdMaster = DS.spectrum.Lambda(resonant_modes);
% % options = struct();
% % options.isauto = true; 
% % options.isdamped = true;
% % options.numDigits = 6;
% % disp('Reduced dynamics on the 2D stable SSM:')
% % y = reduced_dynamics_symbolic(lamdMaster,R0,options)
% DJ{1,i} = lamdMaster;

% hj = matlabFunction(y(2));
% backbonecurve{1,i} = hj;
end    
A = eye(4);
 save('VSP.mat','VT','V','Grid_F');
 save('ASP.mat','AT','A','Grid_F');
 save('DSP.mat','DT','Grid_F');
 save('ICSP.mat','ICTT','Grid_F');

 load('VSP.mat','VT','V','Grid_F');
 load('ASP.mat','AT','A','Grid_F');
 load('DSP.mat','DT','Grid_F');
 load('ICSP.mat','ICTT','Grid_F');


VarA = [];
AarA = [];
DarA = [];
IarA = [];
for ind = 1:max(size(Grid_F))
VarA = [VarA,VT{1,ind}(:)];
AarA = [AarA,full(AT{1,ind}(:))];
DarA = [DarA,DT{1,ind}(:)];
IarA = [IarA,ICTT{1,ind}(:)];
end    


ValphaR = griddedInterpolant(Grid_F(:),real(VarA).','spline');
ValphaI = griddedInterpolant(Grid_F(:),imag(VarA).','spline');
Valpha = @(F,V) reshape((ValphaR(F) + 1i*ValphaI(F)).',max(size(V(:,1))),max(size(V(1,:))));

Aalpha = griddedInterpolant(Grid_F(:),AarA.','spline');
Aalph = @(F,A) reshape(Aalpha(F).',max(size(A(:,1))),max(size(A(1,:))));

Dalpha = griddedInterpolant(Grid_F(:),DarA.','spline');
Dalph = @(F) Dalpha(F);

Ialpha = griddedInterpolant(Grid_F(:),IarA.','spline');
Ialph = @(F) Ialpha(F);

%Order 0
xi_0 = @(alpha) [Ialpha(Force_Lorenz(alpha)).';0;0];
Load_Check = Ialph(F_a_new).';
[dXdF, Xtrunc,tkr1] = ftd(Load_Check, tk.');
IDalpha = griddedInterpolant(tkr1(:),dXdF.','spline');

% [dXdFPa, XtruncPa,F_truncPa] = ftd(F_a_new.', tk.');
%  
% DForce_Lorenz = griddedInterpolant(F_truncPa(:),dXdFPa(:),'spline');

xi_1 = @(alpha,A) (inv(Aalph(Force_Lorenz(alpha),A))*[IDalpha(alpha).';0;0]);
XI_1 = [];
for indf = 1:max(size(tkr1))
    vec1=xi_1(tkr1(indf),A);
    XI_1 = [XI_1,vec1];
end    

% % Order 2 
% 
[dXdF, Xtrunc,tkr2] = ftd(XI_1, tkr1);
xi_1D = griddedInterpolant(tkr2(:),dXdF.','spline');
wrap = @(a,x)a(x);
xi_21 = @(alpha,A) wrap(xi_1(alpha,A),2);
xi_20 = @(alpha) wrap(xi_0(alpha),2);
xi_41 = @(alpha,A) wrap(xi_1(alpha,A),4);

xi_2 = @(alpha,A) (inv(Aalph(Force_Lorenz(alpha),A))*(xi_1D(alpha).'-[0;0;0;(m+Mf)/(m*Mf)*(-4*g/scale*m*3*xi_20(alpha)*xi_21(alpha,A)^2-16*m*a^4/scale^2*xi_20(alpha)*xi_41(alpha,A)^2 +...
   sign*(-4*g/scale*m*3*a*xi_21(alpha,A)^2 -16*m*a^5/scale^2*xi_41(alpha,A)^2) )]));

XI_2 = [];
for indf = 1:max(size(tkr2))
    vec1=xi_2(tkr2(indf),A);
    XI_2 = [XI_2,vec1];
end 

% Order 3
[dXdF, Xtrunc,tkr3] = ftd(XI_2, tkr2);
xi_2D = griddedInterpolant(tkr3(:),dXdF.','spline');
xi_22 = @(alpha,A) wrap(xi_2(alpha,A),2);
xi_42 = @(alpha,A) wrap(xi_2(alpha,A),4);

xi_3 = @(alpha,A) (inv(Aalph(Force_Lorenz(alpha),A))*(xi_2D(alpha).'-[0;0;0;(Mf+m)/(Mf*m)*...
    (-4*g/scale*m*(xi_21(alpha,A)^3 + 6 * xi_20(alpha)*xi_21(alpha,A)*xi_22(alpha,A)) - ...
    16*m*a^4/scale^2*(xi_21(alpha,A)*xi_41(alpha,A)^2 + 2*xi_20(alpha)*xi_41(alpha,A)*xi_42(alpha,A) + ...
    sign*(-4*g/scale*m*3*a*2*xi_21(alpha,A)*xi_22(alpha,A) - 16*m*a^5/scale^2*2*xi_41(alpha,A)*xi_42(alpha,A) ) ))]));

XI_3A = [];
check_2D = [];
for indf = 1:max(size(tkr3))
    vec1=xi_3(tkr3(indf),A);
    check_2D = [check_2D,xi_1D(tkr3(indf)).'];
    XI_3A = [XI_3A,vec1];
end 

XI_2A = [];
for indf = 1:max(size(tkr3))
    vec1=xi_2(tkr3(indf),A);
    XI_2A = [XI_2A,vec1];
end 

XI_1A = [];
for indf = 1:max(size(tkr3))
    vec1=xi_1(tkr3(indf),A);
    XI_1A = [XI_1A,vec1];
end 

XI_0A = [];
for indf = 1:max(size(tkr3))
    vec1=xi_0(tkr3(indf));
    XI_0A = [XI_0A,vec1];
end 

%
XI_0Sp = griddedInterpolant(tkr3(:),XI_0A.','spline');
XI_1Sp = griddedInterpolant(tkr3(:),XI_1A.','spline');
XI_2Sp = griddedInterpolant(tkr3(:),XI_2A.','spline');

porder = 5;
framelen = 101;
sgf1 = sgolayfilt(XI_3A(1,:),porder,framelen);
sgf2 = sgolayfilt(XI_3A(2,:),porder,framelen);
sgf3 = sgolayfilt(XI_3A(3,:),porder,framelen);
sgf4 = sgolayfilt(XI_3A(4,:),porder,framelen);
% sgf5 = sgolayfilt(XI_3A(5,:),porder,framelen);
% sgf6 = sgolayfilt(XI_3A(6,:),porder,framelen);
XI_3f = [sgf1;sgf2;sgf3;sgf4];%sgf5;sgf6];

XI_3Sp = griddedInterpolant(tkr3(:),XI_3f.','spline');
save('Solution_stablep.mat','XI_0Sp','XI_1Sp','XI_2Sp','XI_3Sp');

dias = [];
for ind = 1:1000
    H = VT{1,ind};
    dias = [dias,H(:,4)];
end    
figure 
plot(Grid_F,real(dias(4,:)),'-','LineWidth',3)


dias = [];
for ind = 1:1000
    H = DT{1,ind};
    dias = [dias,H(:)];
end    
figure 
plot(Grid_F,real(dias(1,:)),'-','LineWidth',3)