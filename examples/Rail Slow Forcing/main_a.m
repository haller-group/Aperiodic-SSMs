% [M,C,K,fnl,f_0,outdof,PlotFieldonDefMesh] = build_model(nElements);
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
K = -[-k,k*m/(m+Mf);k*m/(m+Mf),4*a^2*m*g/scale-k*m^2/(m+Mf)^2];
C = -[-c,c*m/(m+Mf);c*m/(m+Mf),-cf-c*m^2/(m+Mf)^2];
M = [m+Mf,0;0,m*Mf/(m+Mf)];



n = 2;
F3 = sptensor([n,n,n,n]);
F2 = sptensor([n,n,n]);
F3(2,2,2,2) = 4*g*m/scale;
%F3(2,2,4,4) = 16*m*a^4;

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
Drot = diag(Ds);
Df = [Drot(4);Drot(1);Drot(2);Drot(3)];
Vs = [Vs(:,4),Vs(:,1),Vs(:,2),Vs(:,3)];

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


% VC = reorient_complex(V,V0);

[Vcheck,Dcheck] = eig(full(AT{1,i}));

dsa = diag(Dcheck);
[d,ind]=sort(real(dsa));
Ds=Dcheck(ind,ind);         %sorted eigenvalue matrix
Vs=Vcheck(:,ind); 
Drot = diag(Ds);
Df = [Drot(4);Drot(1);Drot(2);Drot(3)];
Vs = [Vs(:,4),Vs(:,1),Vs(:,2),Vs(:,3)];

 VC = reorient_complex(Vs,V0);
 inv(VC)
% V0 = VC;
% [V,DA] = eigs(AT{1,i});
% dsa = diag(DA);
% [d,ind]=sort(real(-dsa));
% Ds=DA(ind,ind)          %sorted eigenvalue matrix
% Vs=V(:,ind);
VT{1,i} = VC;
DT{1,i} = Df;

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

xi_2 = @(alpha,A) (inv(Aalph(Force_Lorenz(alpha),A))*(xi_1D(alpha).'-[0;0;0;(m+Mf)/(m*Mf)*(-4*g*m*3*xi_20(alpha)*xi_21(alpha,A)^2/scale-16*m*a^4*xi_20(alpha)*xi_41(alpha,A)^2/scale^2)]));

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
    (-4*g*m/scale*(xi_21(alpha,A)^3 + 6 * xi_20(alpha)*xi_21(alpha,A)*xi_22(alpha,A)) - ...
    16*m*a^4/scale^2*(xi_21(alpha,A)*xi_41(alpha,A)^2 + 2*xi_20(alpha)*xi_41(alpha,A)*xi_42(alpha,A) ))]));

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
XI_0U = griddedInterpolant(tkr3(:),XI_0A.','spline');
XI_1U = griddedInterpolant(tkr3(:),XI_1A.','spline');
XI_2U = griddedInterpolant(tkr3(:),XI_2A.','spline');

porder = 5;
framelen = 101;
sgf1 = sgolayfilt(XI_3A(1,:),porder,framelen);
sgf2 = sgolayfilt(XI_3A(2,:),porder,framelen);
sgf3 = sgolayfilt(XI_3A(3,:),porder,framelen);
sgf4 = sgolayfilt(XI_3A(4,:),porder,framelen);
% sgf5 = sgolayfilt(XI_3A(5,:),porder,framelen);
% sgf6 = sgolayfilt(XI_3A(6,:),porder,framelen);
XI_3f = [sgf1;sgf2;sgf3;sgf4]%;sgf5;sgf6];

XI_3U = griddedInterpolant(tkr3(:),XI_3f.','spline');
save('Solution_Unstable.mat','XI_0U','XI_1U','XI_2U','XI_3U');



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
plot(Grid_F,real(dias(2,:)),'-','LineWidth',3)

% alphaT = linspace(0,6,10000);
% XIU =[];
% XISP=[];
% XISM=[];
% epsilon = 0.01;
% load('Solution_Unstable.mat','XI_0U','XI_1U','XI_2U','XI_3U');
% load('Solution_stablem.mat','XI_0Sm','XI_1Sm','XI_2Sm','XI_3Sm');
% load('Solution_stablep.mat','XI_0Sp','XI_1Sp','XI_2Sp','XI_3Sp');
% 
% for ind = 1:10000
%     alpha=alphaT(ind);
% XI_Unstable = XI_0U(alpha)+ epsilon*XI_1U(alpha)+epsilon^2*XI_2U(alpha)+epsilon^3*XI_3U(alpha);
% XIU =[XIU;XI_Unstable];
% XI_Stablep = XI_0Sp(alpha)+ epsilon*XI_1Sp(alpha)+epsilon^2*XI_2Sp(alpha)+epsilon^3*XI_3Sp(alpha)+ [a*m/(m+Mf),a,0,0];
% XISP =[XISP;XI_Stablep];
% XI_Stablem = XI_0Sm(alpha)+ epsilon*XI_1Sm(alpha)+epsilon^2*XI_2Sm(alpha)+epsilon^3*XI_3Sm(alpha)- [a*m/(m+Mf),a,0,0];
% XISM =[XISM;XI_Stablem];
% end
% 
% plot(alphaT,XIU(:,2),'color','red')
% hold on 
% plot(alphaT,XISP(:,2),'color','blue')
% hold on 
% plot(alphaT,XISM(:,2),'color','blue')
% ctspan = alphaT/(epsilon);
% GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon);
% tic
% IC = (XIU(1,:) + 0.1).';
% [tSP,Full_Traj] = ode45(GH,ctspan,IC);
% toc
% 
% plot(alphaT,XIU(:,2),'color','red')
% hold on 
% plot(alphaT,XISP(:,2),'color','blue')
% hold on 
% plot(alphaT,XISM(:,2),'color','blue')
% hold on 
% plot(alphaT,Full_Traj(:,2),'color','black')





   

% SSM construction 
Coeff_Iterate_1 = [];
for indf = 1:max(size(tkr2))
dia = Dalpha(Force_Lorenz(tkr2(indf)));
l1 = dia(1);
l2 = dia(2);
l3 = dia(3);
alpha = tkr2(indf);
xi20 = xi_20(alpha);

xi21 = xi_21(alpha,A);
xi41 = xi_41(alpha,A);

xi22 = xi_22(alpha,A);
xi42 = xi_42(alpha,A);


P = Valpha(Force_Lorenz(tkr2(indf)),V);
PI = inv(P);


h101 = (-8*(3*g/Sc*m*xi20*xi21*P(2, 1)*PI(3, 4) + 3*g/Sc*Mf*xi20*xi21*P(2, 1)*PI(3, 4) + 4*1/Sc^2*a^4*m*xi20*xi41*P(4, 1)*PI(3, 4) + 4*1/Sc^2*a^4*Mf*xi20*xi41*P(4, 1)*PI(3, 4)))/((l1 - l3)*Mf);
hc101 = conj(h101);
h200 = (-4*(3*g/Sc*m*xi20*P(2, 1)^2*PI(3, 4) + 3*g/Sc*Mf*xi20*P(2, 1)^2*PI(3, 4) + 4*1/Sc^2*a^4*m*xi20*P(4, 1)^2*PI(3, 4) + 4*1/Sc^2*a^4*Mf*xi20*P(4, 1)^2*PI(3, 4)))/((2*l1 - l3)*Mf);
hc200 = conj(h200);

h011 = (-8*(3*g/Sc*m*xi20*xi21*P(2, 2)*PI(3, 4) + 3*g/Sc*Mf*xi20*xi21*P(2, 2)*PI(3, 4) + 4*1/Sc^2*a^4*m*xi20*xi41*P(4, 2)*PI(3, 4) + 4*1/Sc^2*a^4*Mf*xi20*xi41*P(4, 2)*PI(3, 4)))/((l2 - l3)*Mf);
hc011 = conj(h011);

h110 = (-8*(3*g/Sc*m*xi20*P(2, 1)*P(2, 2)*PI(3, 4) + 3*g/Sc*Mf*xi20*P(2, 1)*P(2, 2)*PI(3, 4) + 4*1/Sc^2*a^4*m*xi20*P(4, 1)*P(4, 2)*PI(3, 4) + 4*1/Sc^2*a^4*Mf*xi20*P(4, 1)*P(4, 2)*PI(3, 4)))/ ...
    ((l1 + l2 - l3)*Mf);
hc110 = conj(h110);

h020 = (-4*(3*g/Sc*m*xi20*P(2, 2)^2*PI(3, 4) + 3*g/Sc*Mf*xi20*P(2, 2)^2*PI(3, 4) + 4*1/Sc^2*a^4*m*xi20*P(4, 2)^2*PI(3, 4) + 4*1/Sc^2*a^4*Mf*xi20*P(4, 2)^2*PI(3, 4)))/((2*l2 - l3)*Mf);
hc020 = conj(h020);


h300 = (4*(6*g/Sc*h200*m*xi20*P(2, 1)^2*PI(1, 4) + 6*g/Sc*h200*Mf*xi20*P(2, 1)^2*PI(1, 4) + 8*1/Sc^2*a^4*h200*m*xi20*P(4, 1)^2*PI(1, 4) + 8*1/Sc^2*a^4*h200*Mf*xi20*P(4, 1)^2*PI(1, 4) + ...
      3*g/Sc*h110*m*xi20*P(2, 1)^2*PI(2, 4) + 3*g/Sc*h110*Mf*xi20*P(2, 1)^2*PI(2, 4) + 4*1/Sc^2*a^4*h110*m*xi20*P(4, 1)^2*PI(2, 4) + 4*1/Sc^2*a^4*h110*Mf*xi20*P(4, 1)^2*PI(2, 4) - ...
      g/Sc*m*P(2, 1)^3*PI(3, 4) - g/Sc*Mf*P(2, 1)^3*PI(3, 4) - 6*g/Sc*h200*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g/Sc*h200*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
      6*g/Sc*hc200*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc200*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 4*1/Sc^2*a^4*m*P(2, 1)*P(4, 1)^2*PI(3, 4) - 4*1/Sc^2*a^4*Mf*P(2, 1)*P(4, 1)^2*PI(3, 4) - ...
      8*1/Sc^2*a^4*h200*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h200*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*hc200*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
      8*1/Sc^2*a^4*hc200*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4)))/((3*l1 - l3)*Mf);

h030 = (4*(3*g/Sc*h110*m*xi20*P(2, 2)^2*PI(1, 4) + 3*g/Sc*h110*Mf*xi20*P(2, 2)^2*PI(1, 4) + 4*1/Sc^2*a^4*h110*m*xi20*P(4, 2)^2*PI(1, 4) + 4*1/Sc^2*a^4*h110*Mf*xi20*P(4, 2)^2*PI(1, 4) + ...
      6*g/Sc*h020*m*xi20*P(2, 2)^2*PI(2, 4) + 6*g/Sc*h020*Mf*xi20*P(2, 2)^2*PI(2, 4) + 8*1/Sc^2*a^4*h020*m*xi20*P(4, 2)^2*PI(2, 4) + 8*1/Sc^2*a^4*h020*Mf*xi20*P(4, 2)^2*PI(2, 4) - ...
      g/Sc*m*P(2, 2)^3*PI(3, 4) - g/Sc*Mf*P(2, 2)^3*PI(3, 4) - 6*g/Sc*h020*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g/Sc*h020*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - ...
      6*g/Sc*hc020*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc020*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 4*1/Sc^2*a^4*m*P(2, 2)*P(4, 2)^2*PI(3, 4) - 4*1/Sc^2*a^4*Mf*P(2, 2)*P(4, 2)^2*PI(3, 4) - ...
      8*1/Sc^2*a^4*h020*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h020*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*hc020*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - ...
      8*1/Sc^2*a^4*hc020*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((3*l2 - l3)*Mf);

h120= (4*(6*g/Sc*h110*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 6*g/Sc*h110*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 6*g/Sc*h200*m*xi20*P(2, 2)^2*PI(1, 4) + 6*g/Sc*h200*Mf*xi20*P(2, 2)^2*PI(1, 4) + ...
      8*1/Sc^2*a^4*h110*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 8*1/Sc^2*a^4*h110*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 8*1/Sc^2*a^4*h200*m*xi20*P(4, 2)^2*PI(1, 4) + 8*1/Sc^2*a^4*h200*Mf*xi20*P(4, 2)^2*PI(1, 4) + ...
      12*g/Sc*h020*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 12*g/Sc*h020*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 3*g/Sc*h110*m*xi20*P(2, 2)^2*PI(2, 4) + 3*g/Sc*h110*Mf*xi20*P(2, 2)^2*PI(2, 4) + ...
      16*1/Sc^2*a^4*h020*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 16*1/Sc^2*a^4*h020*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 4*1/Sc^2*a^4*h110*m*xi20*P(4, 2)^2*PI(2, 4) + 4*1/Sc^2*a^4*h110*Mf*xi20*P(4, 2)^2*PI(2, 4) - ...
      3*g/Sc*m*P(2, 1)*P(2, 2)^2*PI(3, 4) - 3*g/Sc*Mf*P(2, 1)*P(2, 2)^2*PI(3, 4) - 6*g/Sc*h020*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g/Sc*h020*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
      6*g/Sc*h110*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g/Sc*h110*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g/Sc*hc020*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
      6*g/Sc*hc020*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc110*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc110*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
      8*1/Sc^2*a^4*m*P(2, 2)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*1/Sc^2*a^4*Mf*P(2, 2)*P(4, 1)*P(4, 2)*PI(3, 4) - 4*1/Sc^2*a^4*m*P(2, 1)*P(4, 2)^2*PI(3, 4) - 4*1/Sc^2*a^4*Mf*P(2, 1)*P(4, 2)^2*PI(3, 4) - ...
      8*1/Sc^2*a^4*h020*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h020*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h110*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
      8*1/Sc^2*a^4*h110*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*hc020*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*1/Sc^2*a^4*hc020*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
      8*1/Sc^2*a^4*hc110*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 8*1/Sc^2*a^4*hc110*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((l1 + 2*l2 - l3)*Mf);


h210 = (4*(3*g/Sc*h110*m*xi20*P(2, 1)^2*PI(1, 4) + 3*g/Sc*h110*Mf*xi20*P(2, 1)^2*PI(1, 4) + 12*g/Sc*h200*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 12*g/Sc*h200*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + ...
      4*1/Sc^2*a^4*h110*m*xi20*P(4, 1)^2*PI(1, 4) + 4*1/Sc^2*a^4*h110*Mf*xi20*P(4, 1)^2*PI(1, 4) + 16*1/Sc^2*a^4*h200*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 16*1/Sc^2*a^4*h200*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + ...
      6*g/Sc*h020*m*xi20*P(2, 1)^2*PI(2, 4) + 6*g/Sc*h020*Mf*xi20*P(2, 1)^2*PI(2, 4) + 6*g/Sc*h110*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 6*g/Sc*h110*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + ...
      8*1/Sc^2*a^4*h020*m*xi20*P(4, 1)^2*PI(2, 4) + 8*1/Sc^2*a^4*h020*Mf*xi20*P(4, 1)^2*PI(2, 4) + 8*1/Sc^2*a^4*h110*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 8*1/Sc^2*a^4*h110*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) - ...
      3*g/Sc*m*P(2, 1)^2*P(2, 2)*PI(3, 4) - 3*g/Sc*Mf*P(2, 1)^2*P(2, 2)*PI(3, 4) - 6*g/Sc*h110*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g/Sc*h110*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
      6*g/Sc*h200*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g/Sc*h200*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g/Sc*hc110*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
      6*g/Sc*hc110*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc200*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g/Sc*hc200*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 4*1/Sc^2*a^4*m*P(2, 2)*P(4, 1)^2*PI(3, 4) - ...
      4*1/Sc^2*a^4*Mf*P(2, 2)*P(4, 1)^2*PI(3, 4) - 8*1/Sc^2*a^4*m*P(2, 1)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*1/Sc^2*a^4*Mf*P(2, 1)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*1/Sc^2*a^4*h110*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
      8*1/Sc^2*a^4*h110*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h200*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*1/Sc^2*a^4*h200*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
      8*1/Sc^2*a^4*hc110*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*1/Sc^2*a^4*hc110*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*1/Sc^2*a^4*hc200*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - ...
      8*1/Sc^2*a^4*hc200*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((2*l1 + l2 - l3)*Mf);



% h101 = (-8*(3*g*m*xi20*xi21*P(2, 1)*PI(3, 4) + 3*g*Mf*xi20*xi21*P(2, 1)*PI(3, 4) + 4*a^4*m*xi20*xi41*P(4, 1)*PI(3, 4) + 4*a^4*Mf*xi20*xi41*P(4, 1)*PI(3, 4)))/((l1 - l3)*Mf);
% hc101 = conj(h101);
% 
% h200 = (-4*(3*g*m*xi20*P(2, 1)^2*PI(3, 4) + 3*g*Mf*xi20*P(2, 1)^2*PI(3, 4) + 4*a^4*m*xi20*P(4, 1)^2*PI(3, 4) + 4*a^4*Mf*xi20*P(4, 1)^2*PI(3, 4)))/((2*l1 - l3)*Mf);
% hc200 = conj(h200);
% 
% h011 = (-8*(3*g*m*xi20*xi21*P(2, 2)*PI(3, 4) + 3*g*Mf*xi20*xi21*P(2, 2)*PI(3, 4) + 4*a^4*m*xi20*xi41*P(4, 2)*PI(3, 4) + 4*a^4*Mf*xi20*xi41*P(4, 2)*PI(3, 4)))/((l2 - l3)*Mf);
% hc011 = conj(h011);
% 
% h110 = (-8*(3*g*m*xi20*P(2, 1)*P(2, 2)*PI(3, 4) + 3*g*Mf*xi20*P(2, 1)*P(2, 2)*PI(3, 4) + 4*a^4*m*xi20*P(4, 1)*P(4, 2)*PI(3, 4) + 4*a^4*Mf*xi20*P(4, 1)*P(4, 2)*PI(3, 4)))/ ...
%     ((l1 + l2 - l3)*Mf);
% hc110 = conj(h110);
% 
% h020 = (-4*(3*g*m*xi20*P(2, 2)^2*PI(3, 4) + 3*g*Mf*xi20*P(2, 2)^2*PI(3, 4) + 4*a^4*m*xi20*P(4, 2)^2*PI(3, 4) + 4*a^4*Mf*xi20*P(4, 2)^2*PI(3, 4)))/((2*l2 - l3)*Mf);
% hc020 = conj(h020);
% 
% 
% h300 = (4*(6*g*h200*m*xi20*P(2, 1)^2*PI(1, 4) + 6*g*h200*Mf*xi20*P(2, 1)^2*PI(1, 4) + 8*a^4*h200*m*xi20*P(4, 1)^2*PI(1, 4) + 8*a^4*h200*Mf*xi20*P(4, 1)^2*PI(1, 4) + ...
%       3*g*h110*m*xi20*P(2, 1)^2*PI(2, 4) + 3*g*h110*Mf*xi20*P(2, 1)^2*PI(2, 4) + 4*a^4*h110*m*xi20*P(4, 1)^2*PI(2, 4) + 4*a^4*h110*Mf*xi20*P(4, 1)^2*PI(2, 4) - ...
%       g*m*P(2, 1)^3*PI(3, 4) - g*Mf*P(2, 1)^3*PI(3, 4) - 6*g*h200*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g*h200*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
%       6*g*hc200*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g*hc200*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 4*a^4*m*P(2, 1)*P(4, 1)^2*PI(3, 4) - 4*a^4*Mf*P(2, 1)*P(4, 1)^2*PI(3, 4) - ...
%       8*a^4*h200*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*a^4*h200*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*a^4*hc200*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
%       8*a^4*hc200*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4)))/((3*l1 - l3)*Mf);
% 
% h030 = (4*(3*g*h110*m*xi20*P(2, 2)^2*PI(1, 4) + 3*g*h110*Mf*xi20*P(2, 2)^2*PI(1, 4) + 4*a^4*h110*m*xi20*P(4, 2)^2*PI(1, 4) + 4*a^4*h110*Mf*xi20*P(4, 2)^2*PI(1, 4) + ...
%       6*g*h020*m*xi20*P(2, 2)^2*PI(2, 4) + 6*g*h020*Mf*xi20*P(2, 2)^2*PI(2, 4) + 8*a^4*h020*m*xi20*P(4, 2)^2*PI(2, 4) + 8*a^4*h020*Mf*xi20*P(4, 2)^2*PI(2, 4) - ...
%       g*m*P(2, 2)^3*PI(3, 4) - g*Mf*P(2, 2)^3*PI(3, 4) - 6*g*h020*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g*h020*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - ...
%       6*g*hc020*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g*hc020*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 4*a^4*m*P(2, 2)*P(4, 2)^2*PI(3, 4) - 4*a^4*Mf*P(2, 2)*P(4, 2)^2*PI(3, 4) - ...
%       8*a^4*h020*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*a^4*h020*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*a^4*hc020*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - ...
%       8*a^4*hc020*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((3*l2 - l3)*Mf);
% 
% h120= (4*(6*g*h110*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 6*g*h110*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 6*g*h200*m*xi20*P(2, 2)^2*PI(1, 4) + 6*g*h200*Mf*xi20*P(2, 2)^2*PI(1, 4) + ...
%       8*a^4*h110*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 8*a^4*h110*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 8*a^4*h200*m*xi20*P(4, 2)^2*PI(1, 4) + 8*a^4*h200*Mf*xi20*P(4, 2)^2*PI(1, 4) + ...
%       12*g*h020*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 12*g*h020*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 3*g*h110*m*xi20*P(2, 2)^2*PI(2, 4) + 3*g*h110*Mf*xi20*P(2, 2)^2*PI(2, 4) + ...
%       16*a^4*h020*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 16*a^4*h020*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 4*a^4*h110*m*xi20*P(4, 2)^2*PI(2, 4) + 4*a^4*h110*Mf*xi20*P(4, 2)^2*PI(2, 4) - ...
%       3*g*m*P(2, 1)*P(2, 2)^2*PI(3, 4) - 3*g*Mf*P(2, 1)*P(2, 2)^2*PI(3, 4) - 6*g*h020*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g*h020*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
%       6*g*h110*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g*h110*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g*hc020*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
%       6*g*hc020*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g*hc110*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g*hc110*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
%       8*a^4*m*P(2, 2)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*a^4*Mf*P(2, 2)*P(4, 1)*P(4, 2)*PI(3, 4) - 4*a^4*m*P(2, 1)*P(4, 2)^2*PI(3, 4) - 4*a^4*Mf*P(2, 1)*P(4, 2)^2*PI(3, 4) - ...
%       8*a^4*h020*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*a^4*h020*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*a^4*h110*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
%       8*a^4*h110*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*a^4*hc020*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*a^4*hc020*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
%       8*a^4*hc110*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 8*a^4*hc110*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((l1 + 2*l2 - l3)*Mf);
% 
% 
% h210 = (4*(3*g*h110*m*xi20*P(2, 1)^2*PI(1, 4) + 3*g*h110*Mf*xi20*P(2, 1)^2*PI(1, 4) + 12*g*h200*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 12*g*h200*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + ...
%       4*a^4*h110*m*xi20*P(4, 1)^2*PI(1, 4) + 4*a^4*h110*Mf*xi20*P(4, 1)^2*PI(1, 4) + 16*a^4*h200*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 16*a^4*h200*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + ...
%       6*g*h020*m*xi20*P(2, 1)^2*PI(2, 4) + 6*g*h020*Mf*xi20*P(2, 1)^2*PI(2, 4) + 6*g*h110*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 6*g*h110*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + ...
%       8*a^4*h020*m*xi20*P(4, 1)^2*PI(2, 4) + 8*a^4*h020*Mf*xi20*P(4, 1)^2*PI(2, 4) + 8*a^4*h110*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + 8*a^4*h110*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) - ...
%       3*g*m*P(2, 1)^2*P(2, 2)*PI(3, 4) - 3*g*Mf*P(2, 1)^2*P(2, 2)*PI(3, 4) - 6*g*h110*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 6*g*h110*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
%       6*g*h200*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g*h200*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 6*g*hc110*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
%       6*g*hc110*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 6*g*hc200*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 6*g*hc200*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 4*a^4*m*P(2, 2)*P(4, 1)^2*PI(3, 4) - ...
%       4*a^4*Mf*P(2, 2)*P(4, 1)^2*PI(3, 4) - 8*a^4*m*P(2, 1)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*a^4*Mf*P(2, 1)*P(4, 1)*P(4, 2)*PI(3, 4) - 8*a^4*h110*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
%       8*a^4*h110*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 8*a^4*h200*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 8*a^4*h200*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
%       8*a^4*hc110*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*a^4*hc110*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 8*a^4*hc200*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - ...
%       8*a^4*hc200*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4)))/((2*l1 + l2 - l3)*Mf);

           

vecC = [h101;h200;h011;h110;h020;h300;h030;h120;h210];

Coeff_Iterate_1 = [Coeff_Iterate_1,vecC];

end
SSM_Coeff_A_1 = griddedInterpolant(tkr2(:),Coeff_Iterate_1.','spline');

[dXdF, Xtrunc,tkr3] = ftd(Coeff_Iterate_1, tkr2);
DSSM_Coeff_A_1 = griddedInterpolant(tkr3(:),dXdF.','spline');

save('Adiabatic_SSM_coeff.mat','SSM_Coeff_A_1')
save('Adiabatic_SSM_coeffD.mat','DSSM_Coeff_A_1')
% 
Coeff_Iterate_2 = [];
for indf = 1:max(size(tkr3))

    dia = Dalpha(Force_Lorenz(tkr3(indf)));
    l1 = dia(1);
    l2 = dia(2);
    l3 = dia(3);
    alpha = tkr3(indf);

    xi20 = xi_20(alpha);

    xi21 = xi_21(alpha,A);
    xi41 = xi_41(alpha,A);

    xi22 = xi_22(alpha,A);
    xi42 = xi_42(alpha,A);

   
    P = Valpha(Force_Lorenz(tkr3(indf)),V);
    PI = inv(P);
    TempCell = num2cell(SSM_Coeff_A_1(alpha));
    [h101,h200,h011,h110,h020,h300,h030,h120,h210]=deal(TempCell{:});
    TempCell = num2cell(DSSM_Coeff_A_1(alpha));
    [hd101,hd200,hd011,hd110,hd020,hd300,hd030,hd120,hd210]=deal(TempCell{:});
    TempCell = num2cell(conj([h101,h200,h011,h110,h020,h300,h030,h120,h210]));
    [hc101,hc200,hc011,hc110,hc020,hc300,hc030,hc120,hc210] = deal(TempCell{:});

    h012 = (-(hd011*Mf) + 24*g/Sc*h101*m*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g/Sc*h101*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 32*1/Sc^2*a^4*h101*m*xi20*xi41*P(4, 2)*PI(1, 4) + ...
     32*1/Sc^2*a^4*h101*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + 24*g/Sc*h011*m*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g/Sc*h011*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 32*1/Sc^2*a^4*h011*m*xi20*xi41*P(4, 2)*PI(2, 4) + ...
     32*1/Sc^2*a^4*h011*Mf*xi20*xi41*P(4, 2)*PI(2, 4) - 12*g/Sc*m*xi21^2*P(2, 2)*PI(3, 4) - 12*g/Sc*Mf*xi21^2*P(2, 2)*PI(3, 4) - 24*g/Sc*m*xi20*xi22*P(2, 2)*PI(3, 4) - ...
     24*g/Sc*Mf*xi20*xi22*P(2, 2)*PI(3, 4) - 16*1/Sc^2*a^4*m*xi41^2*P(2, 2)*PI(3, 4) - 16*1/Sc^2*a^4*Mf*xi41^2*P(2, 2)*PI(3, 4) - 24*g/Sc*h011*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h011*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g/Sc*hc011*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g/Sc*hc011*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi21*xi41*P(4, 2)*PI(3, 4) - ...
     32*1/Sc^2*a^4*Mf*xi21*xi41*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi20*xi42*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi20*xi42*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*h011*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h011*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*hc011*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc011*Mf*xi20*xi41*P(4, 4)*PI(3, 4))/((l2 - l3)*Mf);


h102 = (-(hd101*Mf) + 24*g/Sc*h101*m*xi20*xi21*P(2, 1)*PI(1, 4) + 24*g/Sc*h101*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 32*1/Sc^2*a^4*h101*m*xi20*xi41*P(4, 1)*PI(1, 4) + ...
     32*1/Sc^2*a^4*h101*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 24*g/Sc*h011*m*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g/Sc*h011*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 32*1/Sc^2*a^4*h011*m*xi20*xi41*P(4, 1)*PI(2, 4) + ...
     32*1/Sc^2*a^4*h011*Mf*xi20*xi41*P(4, 1)*PI(2, 4) - 12*g/Sc*m*xi21^2*P(2, 1)*PI(3, 4) - 12*g/Sc*Mf*xi21^2*P(2, 1)*PI(3, 4) - 24*g/Sc*m*xi20*xi22*P(2, 1)*PI(3, 4) - ...
     24*g/Sc*Mf*xi20*xi22*P(2, 1)*PI(3, 4) - 16*1/Sc^2*a^4*m*xi41^2*P(2, 1)*PI(3, 4) - 16*1/Sc^2*a^4*Mf*xi41^2*P(2, 1)*PI(3, 4) - 24*g/Sc*h101*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h101*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g/Sc*hc101*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g/Sc*hc101*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi21*xi41*P(4, 1)*PI(3, 4) - ...
     32*1/Sc^2*a^4*Mf*xi21*xi41*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi20*xi42*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi20*xi42*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*h101*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h101*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*hc101*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc101*Mf*xi20*xi41*P(4, 4)*PI(3, 4))/((l1 - l3)*Mf);
   

h201 = (-(hd200*Mf) + 48*g/Sc*h200*m*xi20*xi21*P(2, 1)*PI(1, 4) + 48*g/Sc*h200*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 12*g/Sc*h101*m*xi20*P(2, 1)^2*PI(1, 4) + ...
     12*g/Sc*h101*Mf*xi20*P(2, 1)^2*PI(1, 4) + 64*1/Sc^2*a^4*h200*m*xi20*xi41*P(4, 1)*PI(1, 4) + 64*1/Sc^2*a^4*h200*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 16*1/Sc^2*a^4*h101*m*xi20*P(4, 1)^2*PI(1, 4) + ... 
     16*1/Sc^2*a^4*h101*Mf*xi20*P(4, 1)^2*PI(1, 4) + 24*g/Sc*h110*m*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g/Sc*h110*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 12*g/Sc*h011*m*xi20*P(2, 1)^2*PI(2, 4) + ...
     12*g/Sc*h011*Mf*xi20*P(2, 1)^2*PI(2, 4) + 32*1/Sc^2*a^4*h110*m*xi20*xi41*P(4, 1)*PI(2, 4) + 32*1/Sc^2*a^4*h110*Mf*xi20*xi41*P(4, 1)*PI(2, 4) + 16*1/Sc^2*a^4*h011*m*xi20*P(4, 1)^2*PI(2, 4) + ...
     16*1/Sc^2*a^4*h011*Mf*xi20*P(4, 1)^2*PI(2, 4) - 12*g/Sc*m*xi21*P(2, 1)^2*PI(3, 4) - 12*g/Sc*Mf*xi21*P(2, 1)^2*PI(3, 4) - 24*g/Sc*h200*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h200*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g/Sc*h101*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g/Sc*h101*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g/Sc*hc200*m*xi20*xi21*P(2, 4)*PI(3, 4) - ...
     24*g/Sc*hc200*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g/Sc*hc101*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g/Sc*hc101*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*m*xi41*P(2, 1)*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi41*P(2, 1)*P(4, 1)*PI(3, 4) - 16*1/Sc^2*a^4*m*xi21*P(4, 1)^2*PI(3, 4) - 16*1/Sc^2*a^4*Mf*xi21*P(4, 1)^2*PI(3, 4) - ...
     32*1/Sc^2*a^4*h200*m*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h200*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h101*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h101*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*hc200*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc200*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*hc101*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc101*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4))/((2*l1 - l3)*Mf);



h021 = (-(hd020*Mf) + 24*g/Sc*h110*m*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g/Sc*h110*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 12*g/Sc*h101*m*xi20*P(2, 2)^2*PI(1, 4) + ...
     12*g/Sc*h101*Mf*xi20*P(2, 2)^2*PI(1, 4) + 32*1/Sc^2*a^4*h110*m*xi20*xi41*P(4, 2)*PI(1, 4) + 32*1/Sc^2*a^4*h110*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + 16*1/Sc^2*a^4*h101*m*xi20*P(4, 2)^2*PI(1, 4) + ...
     16*1/Sc^2*a^4*h101*Mf*xi20*P(4, 2)^2*PI(1, 4) + 48*g/Sc*h020*m*xi20*xi21*P(2, 2)*PI(2, 4) + 48*g/Sc*h020*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 12*g/Sc*h011*m*xi20*P(2, 2)^2*PI(2, 4) + ...
     12*g/Sc*h011*Mf*xi20*P(2, 2)^2*PI(2, 4) + 64*1/Sc^2*a^4*h020*m*xi20*xi41*P(4, 2)*PI(2, 4) + 64*1/Sc^2*a^4*h020*Mf*xi20*xi41*P(4, 2)*PI(2, 4) + 16*1/Sc^2*a^4*h011*m*xi20*P(4, 2)^2*PI(2, 4) + ...
     16*1/Sc^2*a^4*h011*Mf*xi20*P(4, 2)^2*PI(2, 4) - 12*g/Sc*m*xi21*P(2, 2)^2*PI(3, 4) - 12*g/Sc*Mf*xi21*P(2, 2)^2*PI(3, 4) - 24*g/Sc*h020*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h020*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g/Sc*h011*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g/Sc*h011*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g/Sc*hc020*m*xi20*xi21*P(2, 4)*PI(3, 4) - ...
     24*g/Sc*hc020*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g/Sc*hc011*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 24*g/Sc*hc011*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*m*xi41*P(2, 2)*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi41*P(2, 2)*P(4, 2)*PI(3, 4) - 16*1/Sc^2*a^4*m*xi21*P(4, 2)^2*PI(3, 4) - 16*1/Sc^2*a^4*Mf*xi21*P(4, 2)^2*PI(3, 4) - ...
     32*1/Sc^2*a^4*h020*m*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h020*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h011*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h011*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*hc020*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc020*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*hc011*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc011*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4))/((2*l2 - l3)*Mf);


h111 = (-(hd110*Mf) + 24*g/Sc*h110*m*xi20*xi21*P(2, 1)*PI(1, 4) + 24*g/Sc*h110*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 48*g/Sc*h200*m*xi20*xi21*P(2, 2)*PI(1, 4) + ...
     48*g/Sc*h200*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g/Sc*h101*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 24*g/Sc*h101*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 32*1/Sc^2*a^4*h110*m*xi20*xi41*P(4, 1)*PI(1, 4) + ...
     32*1/Sc^2*a^4*h110*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 64*1/Sc^2*a^4*h200*m*xi20*xi41*P(4, 2)*PI(1, 4) + 64*1/Sc^2*a^4*h200*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + ...
     32*1/Sc^2*a^4*h101*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 32*1/Sc^2*a^4*h101*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 48*g/Sc*h020*m*xi20*xi21*P(2, 1)*PI(2, 4) + ...
     48*g/Sc*h020*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g/Sc*h110*m*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g/Sc*h110*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g/Sc*h011*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + ...
     24*g/Sc*h011*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 64*1/Sc^2*a^4*h020*m*xi20*xi41*P(4, 1)*PI(2, 4) + 64*1/Sc^2*a^4*h020*Mf*xi20*xi41*P(4, 1)*PI(2, 4) + ...
     32*1/Sc^2*a^4*h110*m*xi20*xi41*P(4, 2)*PI(2, 4) + 32*1/Sc^2*a^4*h110*Mf*xi20*xi41*P(4, 2)*PI(2, 4) + 32*1/Sc^2*a^4*h011*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + ...
     32*1/Sc^2*a^4*h011*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) - 24*g/Sc*m*xi21*P(2, 1)*P(2, 2)*PI(3, 4) - 24*g/Sc*Mf*xi21*P(2, 1)*P(2, 2)*PI(3, 4) - 24*g/Sc*h110*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h110*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g/Sc*h011*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g/Sc*h011*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
     24*g/Sc*h101*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g/Sc*h101*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g/Sc*hc110*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g/Sc*hc110*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - ...
     24*g/Sc*hc011*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g/Sc*hc011*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g/Sc*hc101*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
     24*g/Sc*hc101*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi41*P(2, 2)*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi41*P(2, 2)*P(4, 1)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi41*P(2, 1)*P(4, 2)*PI(3, 4) - ...
     32*1/Sc^2*a^4*Mf*xi41*P(2, 1)*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*m*xi21*P(4, 1)*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*Mf*xi21*P(4, 1)*P(4, 2)*PI(3, 4) - 32*1/Sc^2*a^4*h110*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h110*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h011*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h011*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
     32*1/Sc^2*a^4*h101*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*h101*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*1/Sc^2*a^4*hc110*m*xi20*xi41*P(4, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*hc110*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc011*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc011*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
     32*1/Sc^2*a^4*hc101*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 32*1/Sc^2*a^4*hc101*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4))/((l1 + l2 - l3)*Mf);



% h012 = (-(hd011*Mf) + 24*g*h101*m*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g*h101*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 32*a^4*h101*m*xi20*xi41*P(4, 2)*PI(1, 4) + ...
%      32*a^4*h101*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + 24*g*h011*m*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g*h011*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 32*a^4*h011*m*xi20*xi41*P(4, 2)*PI(2, 4) + ...
%      32*a^4*h011*Mf*xi20*xi41*P(4, 2)*PI(2, 4) - 12*g*m*xi21^2*P(2, 2)*PI(3, 4) - 12*g*Mf*xi21^2*P(2, 2)*PI(3, 4) - 24*g*m*xi20*xi22*P(2, 2)*PI(3, 4) - ...
%      24*g*Mf*xi20*xi22*P(2, 2)*PI(3, 4) - 16*a^4*m*xi41^2*P(2, 2)*PI(3, 4) - 16*a^4*Mf*xi41^2*P(2, 2)*PI(3, 4) - 24*g*h011*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
%      24*g*h011*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g*hc011*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g*hc011*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 32*a^4*m*xi21*xi41*P(4, 2)*PI(3, 4) - ...
%      32*a^4*Mf*xi21*xi41*P(4, 2)*PI(3, 4) - 32*a^4*m*xi20*xi42*P(4, 2)*PI(3, 4) - 32*a^4*Mf*xi20*xi42*P(4, 2)*PI(3, 4) - 32*a^4*h011*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h011*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*hc011*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*a^4*hc011*Mf*xi20*xi41*P(4, 4)*PI(3, 4))/((l2 - l3)*Mf);
% 
% 
% h102 = (-(hd101*Mf) + 24*g*h101*m*xi20*xi21*P(2, 1)*PI(1, 4) + 24*g*h101*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 32*a^4*h101*m*xi20*xi41*P(4, 1)*PI(1, 4) + ...
%      32*a^4*h101*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 24*g*h011*m*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g*h011*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 32*a^4*h011*m*xi20*xi41*P(4, 1)*PI(2, 4) + ...
%      32*a^4*h011*Mf*xi20*xi41*P(4, 1)*PI(2, 4) - 12*g*m*xi21^2*P(2, 1)*PI(3, 4) - 12*g*Mf*xi21^2*P(2, 1)*PI(3, 4) - 24*g*m*xi20*xi22*P(2, 1)*PI(3, 4) - ...
%      24*g*Mf*xi20*xi22*P(2, 1)*PI(3, 4) - 16*a^4*m*xi41^2*P(2, 1)*PI(3, 4) - 16*a^4*Mf*xi41^2*P(2, 1)*PI(3, 4) - 24*g*h101*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
%      24*g*h101*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g*hc101*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g*hc101*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 32*a^4*m*xi21*xi41*P(4, 1)*PI(3, 4) - ...
%      32*a^4*Mf*xi21*xi41*P(4, 1)*PI(3, 4) - 32*a^4*m*xi20*xi42*P(4, 1)*PI(3, 4) - 32*a^4*Mf*xi20*xi42*P(4, 1)*PI(3, 4) - 32*a^4*h101*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h101*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*hc101*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*a^4*hc101*Mf*xi20*xi41*P(4, 4)*PI(3, 4))/((l1 - l3)*Mf);
%    
% 
% h201 = (-(hd200*Mf) + 48*g*h200*m*xi20*xi21*P(2, 1)*PI(1, 4) + 48*g*h200*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 12*g*h101*m*xi20*P(2, 1)^2*PI(1, 4) + ...
%      12*g*h101*Mf*xi20*P(2, 1)^2*PI(1, 4) + 64*a^4*h200*m*xi20*xi41*P(4, 1)*PI(1, 4) + 64*a^4*h200*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 16*a^4*h101*m*xi20*P(4, 1)^2*PI(1, 4) + ... 
%      16*a^4*h101*Mf*xi20*P(4, 1)^2*PI(1, 4) + 24*g*h110*m*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g*h110*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 12*g*h011*m*xi20*P(2, 1)^2*PI(2, 4) + ...
%      12*g*h011*Mf*xi20*P(2, 1)^2*PI(2, 4) + 32*a^4*h110*m*xi20*xi41*P(4, 1)*PI(2, 4) + 32*a^4*h110*Mf*xi20*xi41*P(4, 1)*PI(2, 4) + 16*a^4*h011*m*xi20*P(4, 1)^2*PI(2, 4) + ...
%      16*a^4*h011*Mf*xi20*P(4, 1)^2*PI(2, 4) - 12*g*m*xi21*P(2, 1)^2*PI(3, 4) - 12*g*Mf*xi21*P(2, 1)^2*PI(3, 4) - 24*g*h200*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
%      24*g*h200*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g*h101*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g*h101*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g*hc200*m*xi20*xi21*P(2, 4)*PI(3, 4) - ...
%      24*g*hc200*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g*hc101*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g*hc101*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - ...
%      32*a^4*m*xi41*P(2, 1)*P(4, 1)*PI(3, 4) - 32*a^4*Mf*xi41*P(2, 1)*P(4, 1)*PI(3, 4) - 16*a^4*m*xi21*P(4, 1)^2*PI(3, 4) - 16*a^4*Mf*xi21*P(4, 1)^2*PI(3, 4) - ...
%      32*a^4*h200*m*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*h200*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*h101*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h101*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 32*a^4*hc200*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*a^4*hc200*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - ...
%      32*a^4*hc101*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 32*a^4*hc101*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4))/((2*l1 - l3)*Mf);
% 
% 
% 
% h021 = (-(hd020*Mf) + 24*g*h110*m*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g*h110*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 12*g*h101*m*xi20*P(2, 2)^2*PI(1, 4) + ...
%      12*g*h101*Mf*xi20*P(2, 2)^2*PI(1, 4) + 32*a^4*h110*m*xi20*xi41*P(4, 2)*PI(1, 4) + 32*a^4*h110*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + 16*a^4*h101*m*xi20*P(4, 2)^2*PI(1, 4) + ...
%      16*a^4*h101*Mf*xi20*P(4, 2)^2*PI(1, 4) + 48*g*h020*m*xi20*xi21*P(2, 2)*PI(2, 4) + 48*g*h020*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 12*g*h011*m*xi20*P(2, 2)^2*PI(2, 4) + ...
%      12*g*h011*Mf*xi20*P(2, 2)^2*PI(2, 4) + 64*a^4*h020*m*xi20*xi41*P(4, 2)*PI(2, 4) + 64*a^4*h020*Mf*xi20*xi41*P(4, 2)*PI(2, 4) + 16*a^4*h011*m*xi20*P(4, 2)^2*PI(2, 4) + ...
%      16*a^4*h011*Mf*xi20*P(4, 2)^2*PI(2, 4) - 12*g*m*xi21*P(2, 2)^2*PI(3, 4) - 12*g*Mf*xi21*P(2, 2)^2*PI(3, 4) - 24*g*h020*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
%      24*g*h020*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g*h011*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g*h011*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g*hc020*m*xi20*xi21*P(2, 4)*PI(3, 4) - ...
%      24*g*hc020*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g*hc011*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 24*g*hc011*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
%      32*a^4*m*xi41*P(2, 2)*P(4, 2)*PI(3, 4) - 32*a^4*Mf*xi41*P(2, 2)*P(4, 2)*PI(3, 4) - 16*a^4*m*xi21*P(4, 2)^2*PI(3, 4) - 16*a^4*Mf*xi21*P(4, 2)^2*PI(3, 4) - ...
%      32*a^4*h020*m*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*h020*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*h011*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h011*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*a^4*hc020*m*xi20*xi41*P(4, 4)*PI(3, 4) - 32*a^4*hc020*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - ...
%      32*a^4*hc011*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 32*a^4*hc011*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4))/((2*l2 - l3)*Mf);
% 
% 
% h111 = (-(hd110*Mf) + 24*g*h110*m*xi20*xi21*P(2, 1)*PI(1, 4) + 24*g*h110*Mf*xi20*xi21*P(2, 1)*PI(1, 4) + 48*g*h200*m*xi20*xi21*P(2, 2)*PI(1, 4) + ...
%      48*g*h200*Mf*xi20*xi21*P(2, 2)*PI(1, 4) + 24*g*h101*m*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 24*g*h101*Mf*xi20*P(2, 1)*P(2, 2)*PI(1, 4) + 32*a^4*h110*m*xi20*xi41*P(4, 1)*PI(1, 4) + ...
%      32*a^4*h110*Mf*xi20*xi41*P(4, 1)*PI(1, 4) + 64*a^4*h200*m*xi20*xi41*P(4, 2)*PI(1, 4) + 64*a^4*h200*Mf*xi20*xi41*P(4, 2)*PI(1, 4) + ...
%      32*a^4*h101*m*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 32*a^4*h101*Mf*xi20*P(4, 1)*P(4, 2)*PI(1, 4) + 48*g*h020*m*xi20*xi21*P(2, 1)*PI(2, 4) + ...
%      48*g*h020*Mf*xi20*xi21*P(2, 1)*PI(2, 4) + 24*g*h110*m*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g*h110*Mf*xi20*xi21*P(2, 2)*PI(2, 4) + 24*g*h011*m*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + ...
%      24*g*h011*Mf*xi20*P(2, 1)*P(2, 2)*PI(2, 4) + 64*a^4*h020*m*xi20*xi41*P(4, 1)*PI(2, 4) + 64*a^4*h020*Mf*xi20*xi41*P(4, 1)*PI(2, 4) + ...
%      32*a^4*h110*m*xi20*xi41*P(4, 2)*PI(2, 4) + 32*a^4*h110*Mf*xi20*xi41*P(4, 2)*PI(2, 4) + 32*a^4*h011*m*xi20*P(4, 1)*P(4, 2)*PI(2, 4) + ...
%      32*a^4*h011*Mf*xi20*P(4, 1)*P(4, 2)*PI(2, 4) - 24*g*m*xi21*P(2, 1)*P(2, 2)*PI(3, 4) - 24*g*Mf*xi21*P(2, 1)*P(2, 2)*PI(3, 4) - 24*g*h110*m*xi20*xi21*P(2, 3)*PI(3, 4) - ...
%      24*g*h110*Mf*xi20*xi21*P(2, 3)*PI(3, 4) - 24*g*h011*m*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - 24*g*h011*Mf*xi20*P(2, 1)*P(2, 3)*PI(3, 4) - ...
%      24*g*h101*m*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g*h101*Mf*xi20*P(2, 2)*P(2, 3)*PI(3, 4) - 24*g*hc110*m*xi20*xi21*P(2, 4)*PI(3, 4) - 24*g*hc110*Mf*xi20*xi21*P(2, 4)*PI(3, 4) - ...
%      24*g*hc011*m*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g*hc011*Mf*xi20*P(2, 1)*P(2, 4)*PI(3, 4) - 24*g*hc101*m*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - ...
%      24*g*hc101*Mf*xi20*P(2, 2)*P(2, 4)*PI(3, 4) - 32*a^4*m*xi41*P(2, 2)*P(4, 1)*PI(3, 4) - 32*a^4*Mf*xi41*P(2, 2)*P(4, 1)*PI(3, 4) - 32*a^4*m*xi41*P(2, 1)*P(4, 2)*PI(3, 4) - ...
%      32*a^4*Mf*xi41*P(2, 1)*P(4, 2)*PI(3, 4) - 32*a^4*m*xi21*P(4, 1)*P(4, 2)*PI(3, 4) - 32*a^4*Mf*xi21*P(4, 1)*P(4, 2)*PI(3, 4) - 32*a^4*h110*m*xi20*xi41*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h110*Mf*xi20*xi41*P(4, 3)*PI(3, 4) - 32*a^4*h011*m*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - 32*a^4*h011*Mf*xi20*P(4, 1)*P(4, 3)*PI(3, 4) - ...
%      32*a^4*h101*m*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*a^4*h101*Mf*xi20*P(4, 2)*P(4, 3)*PI(3, 4) - 32*a^4*hc110*m*xi20*xi41*P(4, 4)*PI(3, 4) - ...
%      32*a^4*hc110*Mf*xi20*xi41*P(4, 4)*PI(3, 4) - 32*a^4*hc011*m*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - 32*a^4*hc011*Mf*xi20*P(4, 1)*P(4, 4)*PI(3, 4) - ...
%      32*a^4*hc101*m*xi20*P(4, 2)*P(4, 4)*PI(3, 4) - 32*a^4*hc101*Mf*xi20*P(4, 2)*P(4, 4)*PI(3, 4))/((l1 + l2 - l3)*Mf);




vecC = [h012;h102;h201;h021;h111];

Coeff_Iterate_2 = [Coeff_Iterate_2,vecC];

end  

SSM_Coeff_A_2 = griddedInterpolant(tkr3(:),Coeff_Iterate_2.','spline');
save('Adiabatic_SSM_coeff.mat','SSM_Coeff_A_1','SSM_Coeff_A_2')

load('Adiabatic_SSM_coeff.mat','SSM_Coeff_A_1','SSM_Coeff_A_2')

epsilon = 0.01;

alphaT = linspace(0,6,10000);
XIU =[];
XISP=[];
XISM=[];

load('Solution_Unstable.mat','XI_0U','XI_1U','XI_2U','XI_3U');
load('Solution_stablem.mat','XI_0Sm','XI_1Sm','XI_2Sm','XI_3Sm');
load('Solution_stablep.mat','XI_0Sp','XI_1Sp','XI_2Sp','XI_3Sp');

for ind = 1:max(size(alphaT))
    alpha=alphaT(ind);
XI_Unstable = XI_0U(alpha)+ epsilon*XI_1U(alpha)+epsilon^2*XI_2U(alpha)+epsilon^3*XI_3U(alpha);
XIU =[XIU;XI_Unstable];
XI_Stablep = XI_0Sp(alpha)+ epsilon*XI_1Sp(alpha)+epsilon^2*XI_2Sp(alpha)+epsilon^3*XI_3Sp(alpha)+ [a*m/(m+Mf),a,0,0];
XISP =[XISP;XI_Stablep];
XI_Stablem = XI_0Sm(alpha)+ epsilon*XI_1Sm(alpha)+epsilon^2*XI_2Sm(alpha)+epsilon^3*XI_3Sm(alpha)- [a*m/(m+Mf),a,0,0];
XISM =[XISM;XI_Stablem];
end

plot(alphaT,XIU(:,2),'color','red')
hold on 
plot(alphaT,XISP(:,2),'color','blue')
hold on 
plot(alphaT,XISM(:,2),'color','blue')

ctspan = alphaT/(epsilon);



ROM=@(t,z) rom_temp_model_adiabatic_rail(t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_20,xi_21,xi_41,xi_22,xi_42,Valpha,V,A,Force_Lorenz,Dalpha,epsilon,a,g,c,cf,Mf,m,k,Sc);

q0 = 0.1 - 0.1*1i;
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,ctspan(1),q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

tic 
[t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
toc

[y,modal,Net_Sol] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_sol,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon,c,k,a,g,m,Mf,cf,scale);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc

% % save('OrderOmega0.1.mat','SP_Traj','y','Net_Sol','tSP');
% 
figure 
hold on 
indexR = 2;
plot(tSP,Full_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(tSP,y(indexR,:),'--','LineWidth',3,'color','red')
hold on 
plot(tSP,XIU(:,indexR),'-','LineWidth',3,'color',[1 0 0 0.3])
hold on 
plot(tSP,XISP(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
hold on 
plot(tSP,XISM(:,indexR),'-','LineWidth',3,'color',[0 0 1 0.3])
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$x \,[$m$]$','Interpreter','latex');




% 
% 
epsilon = 0.01;
 ctspan = linspace(6,20,9000)/epsilon;
dias = [];
epsilon = 0.001;
for indf = 1:max(size(ctspan))
    VF = Valpha(Force_Lorenz(ctspan(indf)*epsilon),V);
    dias = [dias,VF(:,4)];
end    
figure 
plot(ctspan*epsilon,real(dias(4,:)),'-','LineWidth',3)
xlabel('$\alpha = \epsilon t$')
% ylabel('Re[P(1,6)]')


%

% epsilon = 0.001;
%  ctspan = linspace(6,20,90000)/epsilon;
% dias = [];
% epsilon = 0.001;
% for indf = 1:max(size(ctspan))
%     VF = Dalpha(Force_Lorenz(ctspan(indf)*epsilon));
%     dias = [dias,VF.'];
% end    
% figure 
% plot(ctspan*epsilon,real(dias(1,:)),'-','LineWidth',3)
% xlabel('$\alpha = \epsilon t$')
% ylabel('Re[$\lambda_1$]')
%%

% Estimate_ball_cube1 = (alpha*delta/2 - gamma*delta.^3)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-5;
% Estimate_ball_cube2 = (alpha*delta/2 - gamma*delta.^3)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-4;
% Estimate_ball_cube3 = (alpha*delta/2 - gamma*delta.^3)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-3;
% Estimate_ball_cube4 = (alpha*delta/2 - gamma*delta.^3)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-2;
% Estimate_ball_cube5 = (alpha*delta/2 - gamma*delta.^3)./(epsilon*max(abs(F_a_new)));
% 
% Estimate_ball_quad = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% 
% epsilon = 10^-6;
% Estimate_ball_quad1 = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-5;
% Estimate_ball_quad2 = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-4;
% Estimate_ball_quad3 = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-3;
% Estimate_ball_quad4 = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% epsilon = 10^-2;
% Estimate_ball_quad5 = (alpha*delta/2 - gamma*delta.^3-gamma_tilde*delta.^2)./(epsilon*max(abs(F_a_new)));
% 
% 
% figure 
% 
% plot(delta,ones(1,max(size(delta))),'-','LineWidth',3,'color','red');
% hold on 
% plot(delta,Estimate_ball_quad5,'-','LineWidth',1);
% hold on 
% plot(delta,Estimate_ball_quad4,'-','LineWidth',1);
% hold on 
% plot(delta,Estimate_ball_quad3,'-','LineWidth',1);
% hold on 
% plot(delta,Estimate_ball_quad2,'-','LineWidth',1);
% hold on 
% plot(delta,Estimate_ball_quad1,'-','LineWidth',1);
% 
% axis([0,1,-0.5,2])

% epsilon = 10^-2;
% hold on 
% plot(delta,Estimate_ball_quad,'-','LineWidth',3,'color','blue');


% order epsilon solution
% AP = [0,0,1,0;0,0,0,1;-2*k,k,-c,c;k,-2*k,c,-2*c];



% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)));exp(D(5)*(t-ta));conj(exp(D(5)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta)])));
% positive_Lvec = find(tspan<-3999 & tspan>-4000);
% positive_L = positive_Lvec(1)-1;
% 
% 
% tcalc = tspan(positive_L:end);
% N_Max = max(size(tcalc));
% LT = [];
% start_time = [1];
% xi_sol = cell(1,max(size(start_time)));
% 
% tic
% sol_1 = [];
% for i = 1:max(size(tcalc))
% Lo = positive_L - 1 + i ;
% T = tspan(start_time(1):Lo);
% ti = tcalc(i);
% 
% INx = xs(ti,T);
% 
% Q1x = trapz(T,INx,2);
% sol_1  =[sol_1,Q1x];
% end    
% xi_sol{1,1} = sol_1;
% LT = [LT,tspan(start_time(1))];
% toc
%  start_time=positive_L;
% xi_solg=sol_1.'; 
% xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
% xi_full_eval = @(xa) xi_full(xa);
% save('Solution_1.mat','xi_full_eval')
% 
% xi_sol0=xi_sol{1,1};
% xi_sol1 = griddedInterpolant(tcalc(:),(xi_sol0(1,:)).','spline');
% xi_sol2 = @(xq) interp1(tcalc(:),(xi_sol0(2,:)).',xq,'spline');
% xi_sol3 = @(xq) interp1(tcalc(:),(xi_sol0(3,:)).',xq,'spline');
% xi_sol4 = @(xq) interp1(tcalc(:),(xi_sol0(4,:)).',xq,'spline');
% xi_net_sol = @(t) [xi_sol1(t);xi_sol2(t);xi_sol3(t);xi_sol4(t)];
% save('Xi_Sol1.mat','xi_sol1')
% 
% 
% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)));exp(D(5)*(t-ta));conj(exp(D(5)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));-(m1+Mf)/(m1*Mf)*gamma*xi_sol1(ta).^3;1/m1*gamma*xi_sol1(ta).^3;0*ones(1,max(size(ta)))])));
% positive_Lvec = find(tspan<-2999 & tspan>-3000);
% positive_L = positive_Lvec(1)-1;
% % 
% tcalc = tspan(positive_L:end);
% N_Max = max(size(tcalc));
% LT = [];
% 
% xi_sol = cell(1,max(size(start_time)));
% 
% tic
% sol_1 = [];
% for i = 1:max(size(tcalc))
% Lo = positive_L - 1 + i ;
% T = tspan(start_time:Lo);
% ti = tcalc(i);
% 
% INx = xs(ti,T);
% 
% Q1x = trapz(T,INx,2);
% sol_1  =[sol_1,Q1x];
% end    
% xi_sol{1,1} = sol_1;
% LT = [LT,tspan(start_time(1))];
% toc
% 
% xi_solg=sol_1.'; 
% xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
% xi_full_eval3 = @(xa) xi_full(xa);
% save('Solution_3.mat','xi_full_eval3')

%% Calculate time dependent coefficients
% load('Xi_Sol1.mat','xi_sol1')
% load('Solution_3.mat','xi_full_eval3')
% load('Solution_1.mat','xi_full_eval')
% l1 = D(1);
% l1c = conj(D(1));
% l2 = D(3);
% l2c = conj(D(3));
% l3 = D(5);
% l3c = conj(D(5));
% h030 = -(0.013813008195636451 + 0.000996033944032752*1i)*gamma/(l2 - 3*l1c);
% h120 = -(0.04050205667290954 + 0.009257694990541146*1i)*gamma/(l2 - l1-2*l1c);
% h210 = -(0.038622163372940145 + 0.015311760436585583*1i)*gamma/(l2 - 2*l1 - l1c);
% h300 = -(0.011947703447961826 + 0.007003117973441996*1i)*gamma/(l2 - 3*l1);
% 
% f030 = -(0.014143339101206014 - 0.004483904619644111*1i)*gamma/(l3 - 3*l1c);
% f120 = -(0.04398260498843236 - 0.0068399877601822145*1i)*gamma/(l3 - l1-2*l1c);
% f210 = -(0.04451123689935621 - 0.00006902041225595448*1i)*gamma/(l3 - 2*l1 - l1c);
% f300 = -(0.01466786866532256 + 0.0022345179315944852*1i)*gamma/(l3 - 3*l1);
% 
% 
% % load('Excellent_Order_Anchor_solution.mat')
% % xi3_sol1_dummy = @(ta) XI_hyperbolic_solution{1,3}(ta);
% % getc = @(data, cNum) data(:, cNum).';
% % xi3_sol1 = @(xa) getc(xi3_sol1_dummy(xa),1);
% % 
% % ls = D(1);
% % lf = D(3);
% % ap = V(1,1);
% % bp = V(1,3);
% % L1 = VI(1,3);
% % L3 = VI(3,3);
% % L_0 = -4000;
% % tspan_0 = linspace(L_0,300,N_Max);
% % positive_Lvec = find(tspan_0<-2999 & tspan_0>-3000);
% % positive_L = positive_Lvec(1)-1;
% 
% 
% h_coeff=@(t,ta) [-1.*exp((l2-2*l1c)*(t-ta)).*(0.14165458000082587 + 0.021169762212267226*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-2*l1)*(t-ta)).*(0.12873087480450443 + 0.06278965471637281*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-l1-l1c)*(t-ta)).*(0.27356993552351216 + 0.08494825393731548*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-l1c)*(t-ta)).*(0.48134875190217397 + 0.11002354683308162*1i).*gamma.*xi_sol1(ta).^2;-1.*exp((l2-l1)*(t-ta)).*(0.4590070643933819 + 0.1819733954524166*1i).*gamma.*xi_sol1(ta).^2
% ];
% 
% 
% % h021 -1.*exp((l2-2*l1c)*(t-ta)).*(0.14165458000082587 + 0.021169762212267226*1i).*gamma.*xi_sol1(ta)
% % h201 -1.*exp((l2-2*l1)*(t-ta)).*(0.12873087480450443 + 0.06278965471637281*1i).*gamma.*xi_sol1(ta)
% % h111 -1.*exp((l2-l1-l1c)*(t-ta)).*(0.27356993552351216 + 0.08494825393731548*1i).*gamma.*xi_sol1(ta)
% % 
% % h012 -1.*exp((l2-l1c)*(t-ta)).*(0.48134875190217397 + 0.11002354683308162*1i).*gamma.*xi_sol1(ta).^2
% % h102 -1.*exp((l2-l1)*(t-ta)).*(0.4590070643933819 + 0.1819733954524166*1i).*gamma.*xi_sol1(ta).^2
% 
% f_coeff=@(t,ta) [-1.*exp((l3-2*l1c)*(t-ta)).*(0.1493846984852132 - 0.03507901563128705*1i).*gamma.*xi_sol1(ta);-1.*exp((l3-2*l1)*(t-ta)).*(0.1530188850500885 + 0.011469363417641935*1i).*gamma.*xi_sol1(ta);-1.*exp((l3-l1-l1c)*(t-ta)).*(0.3059651596621337 - 0.023887716291801475*1i).*gamma.*xi_sol1(ta);-1.*exp((l3-l1c)*(t-ta)).*(0.5227135053304293 - 0.08129018232281629*1i).*gamma.*xi_sol1(ta).^2;-1.*exp((l3-l1)*(t-ta)).*(0.5289960581546924 - 0.0008202765988769862*1i).*gamma.*xi_sol1(ta).^2];
% 
% % f021 -1.*exp((l3-2*l1c)*(t-ta)).*(0.1493846984852132 - 0.03507901563128705*1i).*gamma.*xi_sol1(ta)
% % f201 -1.*exp((l3-2*l1)*(t-ta)).*(0.1530188850500885 + 0.011469363417641935*1i).*gamma.*xi_sol1(ta)
% % f111 -1.*exp((l3-l1-l1c)*(t-ta)).*(0.3059651596621337 - 0.023887716291801475*1i).*gamma.*xi_sol1(ta)
% % 
% % f012 -1.*exp((l3-l1c)*(t-ta)).*(0.5227135053304293 - 0.08129018232281629*1i).*gamma.*xi_sol1(ta).^2
% % f102 -1.*exp((l3-l1)*(t-ta)).*(0.5289960581546924 - 0.0008202765988769862*1i).*gamma.*xi_sol1(ta).^2
% 
% 
% % tcalc = tspan_0(positive_L:end);
% N_Max = max(size(tcalc));
% sol_1rH = [];
% sol_1iH = [];
% 
% sol_1rF = [];
% sol_1iF = [];
% 
% tic
% for i = 1:max(size(tcalc))
% Lo = positive_L - 1 + i ;
% T = tspan(start_time(1):Lo);
% ti = tcalc(i);
% 
% INxH = h_coeff(ti,T);
% INxrH = real(INxH);
% INxiH = imag(INxH);
% Q1xrH = trapz(T,INxrH,2);
% Q1xiH = trapz(T,INxiH,2);
% sol_1rH  =[sol_1rH,Q1xrH];
% sol_1iH  =[sol_1iH,Q1xiH];
% 
% INxF = f_coeff(ti,T);
% INxrF = real(INxF);
% INxiF = imag(INxH);
% Q1xrF = trapz(T,INxrF,2);
% Q1xiF = trapz(T,INxiF,2);
% sol_1rF  =[sol_1rF,Q1xrF];
% sol_1iF  =[sol_1iF,Q1xiF];
% 
% end    
% toc
% 
% h021r = griddedInterpolant(tcalc(:),(sol_1rH(1,:)).','spline');
% h021i = griddedInterpolant(tcalc(:),(sol_1iH(1,:)).','spline');
% h021 = @(xa) h021r(xa) + 1i*h021i(xa);
% 
% f021r = griddedInterpolant(tcalc(:),(sol_1rF(1,:)).','spline');
% f021i = griddedInterpolant(tcalc(:),(sol_1iF(1,:)).','spline');
% f021 = @(xa) f021r(xa) + 1i*f021i(xa);
% 
% 
% h201r = griddedInterpolant(tcalc(:),(sol_1rH(2,:)).','spline');
% h201i = griddedInterpolant(tcalc(:),(sol_1iH(2,:)).','spline');
% h201 = @(xa) h201r(xa) + 1i*h201i(xa);
% 
% f201r = griddedInterpolant(tcalc(:),(sol_1rF(2,:)).','spline');
% f201i = griddedInterpolant(tcalc(:),(sol_1iF(2,:)).','spline');
% f201 = @(xa) f201r(xa) + 1i*f201i(xa);
% 
% h111r = griddedInterpolant(tcalc(:),(sol_1rH(3,:)).','spline');
% h111i = griddedInterpolant(tcalc(:),(sol_1iH(3,:)).','spline');
% h111 = @(xa) h111r(xa) + 1i*h111i(xa);
% 
% f111r = griddedInterpolant(tcalc(:),(sol_1rF(3,:)).','spline');
% f111i = griddedInterpolant(tcalc(:),(sol_1iF(3,:)).','spline');
% f111 = @(xa) f111r(xa) + 1i*f111i(xa);
% 
% h012r = griddedInterpolant(tcalc(:),(sol_1rH(4,:)).','spline');
% h012i = griddedInterpolant(tcalc(:),(sol_1iH(4,:)).','spline');
% h012 = @(xa) h012r(xa) + 1i*h012i(xa);
% 
% f012r = griddedInterpolant(tcalc(:),(sol_1rF(4,:)).','spline');
% f012i = griddedInterpolant(tcalc(:),(sol_1iF(4,:)).','spline');
% f012 = @(xa) f012r(xa) + 1i*f012i(xa);
% 
% h102r = griddedInterpolant(tcalc(:),(sol_1rH(5,:)).','spline');
% h102i = griddedInterpolant(tcalc(:),(sol_1iH(5,:)).','spline');
% h102 = @(xa) h102r(xa) + 1i*h102i(xa);
% 
% f102r = griddedInterpolant(tcalc(:),(sol_1rF(5,:)).','spline');
% f102i = griddedInterpolant(tcalc(:),(sol_1iF(5,:)).','spline');
% f102 = @(xa) f102r(xa) + 1i*f102i(xa);
% save('Coeff_SSM_NF_Exactv1.mat','h030','h300','h210','h120','h111','h021','h201','h102','h012');
% save('Coeff_SSM_NF_Exactv3.mat','f030','f300','f210','f120','f111','f021','f201','f102','f012');
%%
% load('Coeff_SSM_NF_Exactv1.mat','h030','h300','h210','h120','h111','h021','h201','h102','h012');
% load('Coeff_SSM_NF_Exactv3.mat','f030','f300','f210','f120','f111','f021','f201','f102','f012');
% 
% % Epsilon = [0.1,0.01,0.001,0.0001];
% % NMTET = [];
% % for i = 1:max(size(Epsilon))
% epsilon = 0.01;    
% Order = 3;
% 
% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1);
% ctspan = linspace(0,300,5000);
% 
% q0 = 0.2*exp(1i*0.3);
% y0 =compute_SSM_phy(xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,h111,h021,h201,h102,h012,f030,f300,f210,f120,f111,f021,f201,f102,f012,epsilon);
% 
% IC = [real(q0);imag(q0)];
% 
% tic 
% [t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
% toc
% 
% lf = 0;
% [y,modal] = compute_SSM_phy(xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,h111,h021,h201,h102,h012,f030,f300,f210,f120,f111,f021,f201,f102,f012,epsilon);
% 
% 
% 
% 
% SolutionA = epsilon*xi_full_eval(ctspan.') + epsilon^3*xi_full_eval3(ctspan.');
% FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz);
% 
% tic
% IC = y0;
% [tSP,SP_Traj] = ode45(FullS,ctspan,IC);
% toc
% y = y.';
% figure 
% indexR = 4;
% plot(ctspan,y(:,indexR),'-','LineWidth',3,'color','black')
% hold on 
% plot(ctspan,SP_Traj(:,indexR),'--','LineWidth',3,'color','red')
% hold on 
% plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])
% Rel_Errorz = (sqrt(sum((y - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));
% NMTE = sum(sqrt(sum((y - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));


% NMTET = [NMTET,NMTE];
% end

% figure
% plot(Epsilon,NMTET,'-o','LineWidth',3,'MarkerSize',30,'color','blue')

% figure
% plot(ctspan,Rel_Errorx,'-','LineWidth',3,'color',[0.7 0.7 0.7])
% hold on 
% plot(ctspan,Rel_Error,'-','LineWidth',3,'color',[0 0 0])
% hold on 
% plot(ctspan,Rel_Errorz,'-','LineWidth',3,'color',[0 0.5 0])

% FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz);
% 
% tic
% IC = [0;0;0;0;0;0];
% [tSP,SP_Traj] = ode45(FullS,tcheck,IC);
% toc


