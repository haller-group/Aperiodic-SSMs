
% [M,C,K,fnl,f_0,outdof,PlotFieldonDefMesh] = build_model(nElements);
clear all

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

A = [zeros(3),eye(3);inv(M)*K,inv(M)*C];
% w1 = sqrt(k/m);
% w2 = sqrt(3*k/m);
% ls = (-c/(2*m*w1) + 1i*sqrt(1-c^2/(4*m^2*w1^2)))*sqrt(k/m);
% lf = (-(3*c)/(2*m*w2) + 1i*sqrt(1-(3*c)^2/(4*m^2*w2^2)))*sqrt(3*k/m);

subs2 = [1 1 1 1];
n = 3;
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
set(DS.Options,'Emax',6,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
S = SSM(DS);
resonant_modes = [1 2];
Acheck = V*diag(D)*inv(V);
VI = inv(V);



%Lorenz force 
%%

L = 0;
%L = -5000;
N_Max = 10000;
tspan = linspace(L,500,N_Max);
loren = @(t,y) lorenz_3D(t,y);
IC = [0;0.3;0.5];
[tk,yy] = ode45(loren,tspan,IC);

sigma = 10;
rho = 28;
beta = 8/3;

% F_a = -(sigma*rho*yy(:,1) - sigma*yy(:,1).*yy(:,3) - sigma*yy(:,2) - sigma^2*yy(:,2) + sigma^2*yy(:,1));
F_a_new = yy(:,1)/max(abs(yy(:,1)));


Force_Lorenz = griddedInterpolant(tspan(:),F_a_new(:),'spline');



fig1 = figure
plot(tk,F_a_new,'-','LineWidth',3,'Color','red')
xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$F(t) \,[$N$]$','Interpreter','latex');
axis([0 50 -1 1])


% set(fig1, 'Units' , 'Inches' );
% pos = get(fig1, 'Position' );
% set(fig1, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig1, 'forcing_weak.pdf' , '-dpdf' , '-r300' )


% pos_t = find(tk>0);
% figure
% plot(tk(pos_t),F_a_new(pos_t),'Color','red')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$F(t) \,[$N$]$','Interpreter','latex');
% axis([0 300 -1 1])

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


%%
lambda = [D(1);conj(D(1));D(3);conj(D(3));D(5);conj(D(5))];
positive_L = 1;
tcalc = tspan(positive_L:end);
Ng = 1;
[xG,wG]=lgwt(Ng,0,tcalc(2)-tcalc(1));
tG_new = flip(xG) + tcalc ;
TGnew = tG_new(:);
TGnew = TGnew(1:end-Ng).';
etG = tcalc - tcalc(1);
t_minus_s = (tcalc(2)-tcalc(1))-flip(xG);
etG_minus = flip(t_minus_s) + etG;
etGM= etG_minus(:);
etGM = etGM(1:end-Ng).';
WG = repmat(wG.',1,max(size(tcalc))-1);

FtG  = [sparse(2*n-1,max(size(TGnew)));Force_Lorenz(TGnew.').'];
phiG = VI*FtG;% (rearrange IMP)
GreenPos = exp(lambda.*etGM); 
%               GreenPos = GreenPos(find(GreenPos~=0)); 
Explambda = GreenPos;
eta = zeros(2*size(M,1),max(size(tcalc)));
for j = 1:2*size(M,1)  
    eta_d = conv(WG.*phiG(j,:),Explambda(j,:));
    ETA_M = eta_d(Ng:Ng:max(size(TGnew)));
    eta(j,:) = [0,ETA_M];
end
z = real(V*eta);
xi_solg = z.';
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval = @(xa) xi_full(xa);
xi_sol1 = griddedInterpolant(tcalc(:),(z(1,:)).','spline');
save('Solution_1.mat','xi_full_eval')

% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)));exp(D(5)*(t-ta));conj(exp(D(5)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta)])));
% 
% N_Max = max(size(tcalc));
% LT = [];
%  start_time = [1];
% xi_sol = cell(1,max(size(start_time)));
% 
% tic
% sol_1 = [];
% parfor i = 1:max(size(tcalc))
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
% xi_sol1 = griddedInterpolant(tcalc(:),(z(1,:)).','spline');
% % xi_sol2 = @(xq) interp1(tcalc(:),(xi_sol0(2,:)).',xq,'spline');
% % xi_sol3 = @(xq) interp1(tcalc(:),(xi_sol0(3,:)).',xq,'spline');
% % xi_sol4 = @(xq) interp1(tcalc(:),(xi_sol0(4,:)).',xq,'spline');
% % xi_net_sol = @(t) [xi_sol1(t);xi_sol2(t);xi_sol3(t);xi_sol4(t)];
% save('Xi_Sol1.mat','xi_sol1')


% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)));exp(D(5)*(t-ta));conj(exp(D(5)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));-(m1+Mf)/(m1*Mf)*gamma*xi_sol1(ta).^3;1/m1*gamma*xi_sol1(ta).^3;0*ones(1,max(size(ta)))])));
% % positive_Lvec = find(tspan<-2999 & tspan>-3000);
% % positive_L = positive_Lvec(1)-1;
% positive_L = 1;
% 
% tcalc = tspan(positive_L:end);
% N_Max = max(size(tcalc));
% LT = [];
% 
% xi_sol = cell(1,max(size(1)));
% 
% tic
% sol_1 = [];
% parfor i = 1:max(size(tcalc))
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
% 
% xi_sol0=xi_sol{1,1};
% xi_sol3 = griddedInterpolant(tcalc(:),(xi_sol0(1,:)).','spline');
% save('Xi_Sol3.mat','xi_sol3')

FtG  = [0*ones(1,max(size(TGnew)));0*ones(1,max(size(TGnew)));0*ones(1,max(size(TGnew)));-(m1+Mf)/(m1*Mf)*gamma*xi_sol1(TGnew).^3;1/m1*gamma*xi_sol1(TGnew).^3;0*ones(1,max(size(TGnew)))];
phiG = VI*FtG;% (rearrange IMP)
% GreenPos = exp(lambda.*etGM); 
%               GreenPos = GreenPos(find(GreenPos~=0)); 
% Explambda = GreenPos;
eta = zeros(2*size(M,1),max(size(tcalc)));
for j = 1:2*size(M,1)  
    eta_d = conv(WG.*phiG(j,:),Explambda(j,:));
    ETA_M = eta_d(Ng:Ng:max(size(TGnew)));
    eta(j,:) = [0,ETA_M];
end
z = real(V*eta);
xi_solg = z.';
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval3 = @(xa) xi_full(xa);
xi_sol3 = griddedInterpolant(tcalc(:),(z(1,:)).','spline');
save('Solution_3.mat','xi_full_eval3')
save('Xi_Sol3.mat','xi_sol3')



% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)));exp(D(5)*(t-ta));conj(exp(D(5)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));3*-(m1+Mf)/(m1*Mf)*gamma*xi_sol1(ta).^2.*xi_sol3(ta);3*1/m1*gamma*xi_sol1(ta).^2.*xi_sol3(ta);0*ones(1,max(size(ta)))])));
% % positive_Lvec = find(tspan<-2999 & tspan>-3000);
% % positive_L = positive_Lvec(1)-1;
% positive_L = 1;
% 
% tcalc = tspan(positive_L:end);
% N_Max = max(size(tcalc));
% LT = [];
% 
% xi_sol = cell(1,max(size(start_time)));
% 
% tic
% sol_1 = [];
% parfor i = 1:max(size(tcalc))
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
% xi_full_eval5 = @(xa) xi_full(xa);
% save('Solution_5.mat','xi_full_eval5')


FtG  = [0*ones(1,max(size(TGnew)));0*ones(1,max(size(TGnew)));0*ones(1,max(size(TGnew)));3*(-(m1+Mf)/(m1*Mf))*gamma*xi_sol1(TGnew).^2.*xi_sol3(TGnew);3*1/m1*gamma*xi_sol1(TGnew).^2.*xi_sol3(TGnew);0*ones(1,max(size(TGnew)))];
phiG = VI*FtG;% (rearrange IMP)
% GreenPos = exp(lambda.*etGM); 
%               GreenPos = GreenPos(find(GreenPos~=0)); 
% Explambda = GreenPos;
eta = zeros(2*size(M,1),max(size(tcalc)));
for j = 1:2*size(M,1)  
    eta_d = conv(WG.*phiG(j,:),Explambda(j,:));
    ETA_M = eta_d(Ng:Ng:max(size(TGnew)));
    eta(j,:) = [0,ETA_M];
end
z = real(V*eta);
xi_solg = z.';
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval5 = @(xa) xi_full(xa);
xi_sol5 = griddedInterpolant(tcalc(:),(z(1,:)).','spline');
save('Solution_5.mat','xi_full_eval5')
save('Xi_Sol5.mat','xi_sol5')



%% Calculate time dependent coefficients
load('Xi_Sol1.mat','xi_sol1')
load('Solution_3.mat','xi_full_eval3')
load('Solution_1.mat','xi_full_eval')
load('Solution_5.mat','xi_full_eval5')
load('Xi_Sol3.mat','xi_sol3')
l1 = D(1);
l1c = conj(D(1));
l2 = D(3);
l2c = conj(D(3));
l3 = D(5);
l3c = conj(D(5));
% h030 = -(0.013813008195636451 + 0.000996033944032752*1i)*gamma/(l2 - 3*l1c);
% h120 = -(0.04050205667290954 + 0.009257694990541146*1i)*gamma/(l2 - l1-2*l1c);
% h210 = -(0.038622163372940145 + 0.015311760436585583*1i)*gamma/(l2 - 2*l1 - l1c);
% h300 = -(0.011947703447961826 + 0.007003117973441996*1i)*gamma/(l2 - 3*l1);
% 
% f030 = -(0.014143339101206014 - 0.004483904619644111*1i)*gamma/(l3 - 3*l1c);
% f120 = -(0.04398260498843236 - 0.0068399877601822145*1i)*gamma/(l3 - l1-2*l1c);
% f210 = -(0.04451123689935621 - 0.00006902041225595448*1i)*gamma/(l3 - 2*l1 - l1c);
% f300 = -(0.01466786866532256 + 0.0022345179315944852*1i)*gamma/(l3 - 3*l1);


[h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ] = Auto_Coeffs(gamma,l2,l3,l1c,l1);

% load('Excellent_Order_Anchor_solution.mat')
% xi3_sol1_dummy = @(ta) XI_hyperbolic_solution{1,3}(ta);
% getc = @(data, cNum) data(:, cNum).';
% xi3_sol1 = @(xa) getc(xi3_sol1_dummy(xa),1);
% 
% ls = D(1);
% lf = D(3);
% ap = V(1,1);
% bp = V(1,3);
% L1 = VI(1,3);
% L3 = VI(3,3);
% L_0 = -4000;
% tspan_0 = linspace(L_0,300,N_Max);
% positive_Lvec = find(tspan_0<-2999 & tspan_0>-3000);
% positive_L = positive_Lvec(1)-1;

%%
h_coeff = @(t,ta) H_coeff_1(t,ta,gamma,l1,l2,l3,l1c,h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ,xi_sol1,xi_sol3);
% h_coeff=@(t,ta) [-1.*exp((l2-2*l1c)*(t-ta)).*(0.14165458000082587 + 0.021169762212267226*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-2*l1)*(t-ta)).*(0.12873087480450443 + 0.06278965471637281*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-l1-l1c)*(t-ta)).*(0.27356993552351216 + 0.08494825393731548*1i).*gamma.*xi_sol1(ta);-1.*exp((l2-l1c)*(t-ta)).*(0.48134875190217397 + 0.11002354683308162*1i).*gamma.*xi_sol1(ta).^2;-1.*exp((l2-l1)*(t-ta)).*(0.4590070643933819 + 0.1819733954524166*1i).*gamma.*xi_sol1(ta).^2
% ];
% 
% 
% % h021 -1.*exp((l2-2*l1c)*(t-ta)).*(0.14165458000082587 + 0.021169762212267226*1i).*gamma.*xi_sol1(ta)
% % h201 -1.*exp((l2-2*l1)*(t-ta)).*(0.12873087480450443 + 0.06278965471637281*1i).*gamma.*xi_sol1(ta)
% % h111 -1.*exp((l2-l1-l1c)*(t-ta)).*(0.27356993552351216 +
% % 0.08494825393731548*1i).*gamma.*xi_sol1(ta)   
% % 
% % h012 -1.*exp((l2-l1c)*(t-ta)).*(0.48134875190217397 + 0.11002354683308162*1i).*gamma.*xi_sol1(ta).^2
% % h102 -1.*exp((l2-l1)*(t-ta)).*(0.4590070643933819 + 0.1819733954524166*1i).*gamma.*xi_sol1(ta).^2
% 
% 

f_coeff = @(t,ta) F_coeff_1(t,ta,gamma,l1,l2,l3,l1c,h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ,xi_sol1,xi_sol3);

% % f021 -1.*exp((l3-2*l1c)*(t-ta)).*(0.1493846984852132 - 0.03507901563128705*1i).*gamma.*xi_sol1(ta)
% % f201 -1.*exp((l3-2*l1)*(t-ta)).*(0.1530188850500885 + 0.011469363417641935*1i).*gamma.*xi_sol1(ta)
% % f111 -1.*exp((l3-l1-l1c)*(t-ta)).*(0.3059651596621337 - 0.023887716291801475*1i).*gamma.*xi_sol1(ta)
% % 
% % f012 -1.*exp((l3-l1c)*(t-ta)).*(0.5227135053304293 - 0.08129018232281629*1i).*gamma.*xi_sol1(ta).^2
% % f102 -1.*exp((l3-l1)*(t-ta)).*(0.5289960581546924 - 0.0008202765988769862*1i).*gamma.*xi_sol1(ta).^2
% 
% 

% exp((l2-2.*l1c).*(t-ta))
% exp((l2-2.*l1).*(t-ta))
% exp((l2-l1c-l1).*(t-ta))
% exp((l2-l1c).*(t-ta))
% exp((l2-l1).*(t-ta))

% exp((l3-2.*l1c).*(t-ta)).*
% exp((l3-2.*l1).*(t-ta)).*
% exp((l3-l1c-l1).*(t-ta)).*
% exp((l3-l1c).*(t-ta)).*
% exp((l3-l1).*(t-ta)).*

lambdah = [l2-2*l1c;l2-2*l1;l2-l1c-l1;l2-l1c;l2-l1]; 
GreenPos = exp(lambdah.*etGM); 
Explambdah = GreenPos;

lambdaf = [l3-2.*l1c;l3-2.*l1;l3-l1c-l1;l3-l1c;l3-l1];
GreenPos = exp(lambdaf.*etGM); 
Explambdaf = GreenPos;


phiGh = h_coeff(0,TGnew);% (rearrange IMP)
phiGf = f_coeff(0,TGnew);% (rearrange IMP)
etah = zeros(size(phiGh,1),max(size(tcalc)));
etaf = zeros(size(phiGf,1),max(size(tcalc)));
for j = 1:size(phiGh,1)  
    eta_dh = conv(WG.*phiGh(j,:),Explambdah(j,:));
    eta_df = conv(WG.*phiGf(j,:),Explambdaf(j,:));
    ETA_Mh = eta_dh(Ng:Ng:max(size(TGnew)));
    ETA_Mf = eta_df(Ng:Ng:max(size(TGnew)));
    etah(j,:) = [0,ETA_Mh];
    etaf(j,:) = [0,ETA_Mf];
end
H_Coeff1=griddedInterpolant(tcalc(:),etah.','spline');
F_Coeff1=griddedInterpolant(tcalc(:),etaf.','spline');
save('Coeff_SSM_Exact_O1.mat','H_Coeff1','F_Coeff1');

load('Coeff_SSM_Exact_O1.mat','H_Coeff1','F_Coeff1');


% positive_L = 1;
% tcalc = tspan(1:end);
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
% INxiF = imag(INxF);
% Q1xrF = trapz(T,INxrF,2);
% Q1xiF = trapz(T,INxiF,2);
% sol_1rF  =[sol_1rF,Q1xrF];
% sol_1iF  =[sol_1iF,Q1xiF];
% 
% end    
% toc
% sol_1H = sol_1rH + 1i*sol_1iH;
% sol_1F = sol_1rF + 1i*sol_1iF;
% 
% H_Coeff1=griddedInterpolant(tcalc(:),sol_1H.','spline');
% F_Coeff1=griddedInterpolant(tcalc(:),sol_1F.','spline');
% save('Coeff_SSM_Exact_O1.mat','H_Coeff1','F_Coeff1');
% 
% load('Coeff_SSM_Exact_O1.mat','H_Coeff1','F_Coeff1');

% % h021r = griddedInterpolant(tcalc(:),(sol_1rH(1,:)).','spline');
% % h021i = griddedInterpolant(tcalc(:),(sol_1iH(1,:)).','spline');
% % h021 = @(xa) h021r(xa) + 1i*h021i(xa);
% % 
% % f021r = griddedInterpolant(tcalc(:),(sol_1rF(1,:)).','spline');
% % f021i = griddedInterpolant(tcalc(:),(sol_1iF(1,:)).','spline');
% % f021 = @(xa) f021r(xa) + 1i*f021i(xa);
% % 
% % 
% % h201r = griddedInterpolant(tcalc(:),(sol_1rH(2,:)).','spline');
% % h201i = griddedInterpolant(tcalc(:),(sol_1iH(2,:)).','spline');
% % h201 = @(xa) h201r(xa) + 1i*h201i(xa);
% % 
% % f201r = griddedInterpolant(tcalc(:),(sol_1rF(2,:)).','spline');
% % f201i = griddedInterpolant(tcalc(:),(sol_1iF(2,:)).','spline');
% % f201 = @(xa) f201r(xa) + 1i*f201i(xa);
% % 
% % h111r = griddedInterpolant(tcalc(:),(sol_1rH(3,:)).','spline');
% % h111i = griddedInterpolant(tcalc(:),(sol_1iH(3,:)).','spline');
% % h111 = @(xa) h111r(xa) + 1i*h111i(xa);
% % 
% % f111r = griddedInterpolant(tcalc(:),(sol_1rF(3,:)).','spline');
% % f111i = griddedInterpolant(tcalc(:),(sol_1iF(3,:)).','spline');
% % f111 = @(xa) f111r(xa) + 1i*f111i(xa);
% % 
% % h012r = griddedInterpolant(tcalc(:),(sol_1rH(4,:)).','spline');
% % h012i = griddedInterpolant(tcalc(:),(sol_1iH(4,:)).','spline');
% % h012 = @(xa) h012r(xa) + 1i*h012i(xa);
% % 
% % f012r = griddedInterpolant(tcalc(:),(sol_1rF(4,:)).','spline');
% % f012i = griddedInterpolant(tcalc(:),(sol_1iF(4,:)).','spline');
% % f012 = @(xa) f012r(xa) + 1i*f012i(xa);
% % 
% % h102r = griddedInterpolant(tcalc(:),(sol_1rH(5,:)).','spline');
% % h102i = griddedInterpolant(tcalc(:),(sol_1iH(5,:)).','spline');
% % h102 = @(xa) h102r(xa) + 1i*h102i(xa);
% % 
% % f102r = griddedInterpolant(tcalc(:),(sol_1rF(5,:)).','spline');
% % f102i = griddedInterpolant(tcalc(:),(sol_1iF(5,:)).','spline');
% % f102 = @(xa) f102r(xa) + 1i*f102i(xa);
% 
% 
% 
% 
% 
% % save('Coeff_SSM_NF_Exactv1.mat','h030','h300','h210','h120','h111','h021','h201','h102','h012');
% % save('Coeff_SSM_NF_Exactv3.mat','f030','f300','f210','f120','f111','f021','f201','f102','f012');
% % 
% % load('Coeff_SSM_NF_Exactv1.mat','h030','h300','h210','h120','h111','h021','h201','h102','h012');
% % load('Coeff_SSM_NF_Exactv3.mat','f030','f300','f210','f120','f111','f021','f201','f102','f012');
% 
% 


% exp((l2-4.*l1).*(t-ta)).*
% exp((l2-4.*l1c).*(t-ta)).*
% exp((l2-l1c-3.*l1).*(t-ta)).*

% exp((l2-3.*l1c-l1).*(t-ta)).*

% exp((l2-2.*l1c-2.*l1).*(t-ta)).*
% exp((l2-3.*l1).*(t-ta)).*
% exp((l2-3.*l1c).*(t-ta)).*
% exp((l2-l1c-2.*l1).*(t-ta)).*
% exp((l2-2.*l1c-l1).*(t-ta)).*
% exp((l2-2.*l1).*(t-ta)).*
% exp((l2-2.*l1c).*(t-ta)).*
% exp((l2-l1c-l1).*(t-ta)).*
% exp((l2-l1).*(t-ta)).*
% exp((l2-l1c).*(t-ta)).*

h_coeff_2 = @(t,ta) H_coeff_2(t,ta,gamma,l1,l2,l3,l1c,h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ,H_Coeff1,F_Coeff1,xi_sol1,xi_sol3);

f_coeff_2 = @(t,ta) F_coeff_2(t,ta,gamma,l1,l2,l3,l1c,h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ,H_Coeff1,F_Coeff1,xi_sol1,xi_sol3);

lambdah = [l2-4.*l1;l2-4.*l1c;l2-l1c-3.*l1;l2-3.*l1c-l1;
    l2-2.*l1c-2.*l1;l2-3.*l1;l2-3.*l1c;l2-l1c-2.*l1;l2-2.*l1c-l1;l2-2.*l1;
    l2-2.*l1c;l2-l1c-l1;l2-l1;l2-l1c]; 

GreenPos = exp(lambdah.*etGM); 
Explambdah = GreenPos;

lambdaf = [l3-4.*l1;l3-4.*l1c;l3-l1c-3.*l1;l3-3.*l1c-l1;l3-2.*l1c-2.*l1;l3-3.*l1;l3-3.*l1c;l3-l1c-2.*l1;l3-2.*l1c-l1;l3-2.*l1;l3-2.*l1c;l3-l1c-l1;l3-l1;l3-l1c]; 

GreenPos = exp(lambdaf.*etGM); 
Explambdaf = GreenPos;


phiGh = h_coeff_2(0,TGnew);% (rearrange IMP)
phiGf = f_coeff_2(0,TGnew);% (rearrange IMP)
etah = zeros(size(phiGh,1),max(size(tcalc)));
etaf = zeros(size(phiGf,1),max(size(tcalc)));
for j = 1:size(phiGh,1)  
    eta_dh = conv(WG.*phiGh(j,:),Explambdah(j,:));
    eta_df = conv(WG.*phiGf(j,:),Explambdaf(j,:));
    ETA_Mh = eta_dh(Ng:Ng:max(size(TGnew)));
    ETA_Mf = eta_df(Ng:Ng:max(size(TGnew)));
    etah(j,:) = [0,ETA_Mh];
    etaf(j,:) = [0,ETA_Mf];
end
H_Coeff2=griddedInterpolant(tcalc(:),etah.','spline');
F_Coeff2=griddedInterpolant(tcalc(:),etaf.','spline');
save('Coeff_SSM_Exact_O2.mat','H_Coeff2','F_Coeff2');

load('Coeff_SSM_Exact_O2.mat','H_Coeff2','F_Coeff2');



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
% INxH = h_coeff_2(ti,T);
% INxrH = real(INxH);
% INxiH = imag(INxH);
% Q1xrH = trapz(T,INxrH,2);
% Q1xiH = trapz(T,INxiH,2);
% sol_1rH  =[sol_1rH,Q1xrH];
% sol_1iH  =[sol_1iH,Q1xiH];
% 
% INxF = f_coeff_2(ti,T);
% INxrF = real(INxF);
% INxiF = imag(INxF);
% Q1xrF = trapz(T,INxrF,2);
% Q1xiF = trapz(T,INxiF,2);
% sol_1rF  =[sol_1rF,Q1xrF];
% sol_1iF  =[sol_1iF,Q1xiF];
% 
% end    
% toc
% 
% sol_1H = sol_1rH + 1i*sol_1iH;
% sol_1F = sol_1rF + 1i*sol_1iF;
% 
% H_Coeff2=griddedInterpolant(tcalc(:),sol_1H.','spline');
% F_Coeff2=griddedInterpolant(tcalc(:),sol_1F.','spline');
% save('Coeff_SSM_Exact_O2.mat','H_Coeff2','F_Coeff2');
% 
% load('Coeff_SSM_Exact_O2.mat','H_Coeff2','F_Coeff2');

% Epsilon = [0.1,0.01,0.001,0.0001];
% NMTET = [];
% for i = 1:max(size(Epsilon))

% epsilon = 0;    
% Order = 3;
% 

%%
epsilon = 0.5;

ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);
ROM_ODE_linear =@(t,zk) rom_temp_model_linear(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);

ctspan = linspace(0,500,8000);
%dy = rom_temp_model_linear(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1)

q0 = 0.6*exp(1i*0);
y0 =compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);

q0 = 0.6*exp(1i*0);
y0_l =compute_SSM_linear(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);


IC = [real(q0);imag(q0)];

% IC = [1.2;0];
tic 
[t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
toc

tic 
[t_sol,ROM_sol_linear] = ode45(ROM_ODE_linear,ctspan,IC);
toc

lf = 0;
[y,modal] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);

[y_linear,modal_linear] = compute_SSM_linear(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol_linear(:,1)+1i*ROM_sol_linear(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);



SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);

tic
IC =  y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc
y = y.';
figure 
indexR = 3;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
% plot(ctspan,y_linear(indexR,:),'--','LineWidth',3,'color','blue')
% hold on
plot(ctspan,SolutionA(:,indexR),'-','LineWidth',3,'color',[0 1 0 0.3])
 

Rel_Errorz = (sqrt(sum((y - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2))));
NMTE = sum(sqrt(sum((y - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));

%%
% Epsilon = [0.001,0.005,0.01,0.05,0.1,0.5];
% NMTET = [];
% for ij = 1:max(size(Epsilon))
%  epsilon = Epsilon(ij);
% 
% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1,xi_sol3,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1);
% ctspan = linspace(0,500,8000);
% 
% q0 = 1.2*exp(1i*0);
% y0 =compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,0,l1,l1,gamma,V,inv(V),q0,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);
% 
% IC = [real(q0);imag(q0)];
% 
% % IC = [1.2;0];
% tic 
% [t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
% toc
% 
% lf = 0;
% [y,modal] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);
% 
% 
% 
% 
% SolutionA = epsilon*xi_full_eval(ctspan.') +epsilon^3*xi_full_eval3(ctspan.')+epsilon^5*xi_full_eval5(ctspan.');
% FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);
% 
% tic
% IC =  y0;
% [tSP,SP_Traj] = ode45(FullS,ctspan,IC);
% toc
% y = y.';   
% NMTE = sum(sqrt(sum((y - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));
% NMTET = [NMTET,NMTE];
% end
% fig = figure 
% plot((m1+m2+Mf)*Epsilon,NMTET,'-o','MarkerSize',10,'LineWidth',3,'color','blue')
% xlabel('max$(|F_{ext}|)$ [$\mathbf{N}$]','Interpreter','latex');
% ylabel('NMTE','Interpreter','latex');
% title({'Normalised mean trajectory error' 'for varying max$(|F_{ext}|)$'},'Interpreter','latex')

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'NMTE_weak.pdf' , '-dpdf' , '-r300' )

% title_string = strcat('max$(|F_{ext}|) = ',num2str(epsilon),'$ [$\mathbf{N}$], expansion order $N = 5$');
% title(ax1,title_string,'Interpreter','latex')

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

% 
% h021r = griddedInterpolant(tcalc(:),(sol_1r(2,:)).','spline');
% h021i = griddedInterpolant(tcalc(:),(sol_1i(2,:)).','spline');
% h021 = @(xa) h021r(xa) + 1i*h021i(xa);
% 
% 
% h111r = griddedInterpolant(tcalc(:),(sol_1r(3,:)).','spline');
% h111i = griddedInterpolant(tcalc(:),(sol_1i(3,:)).','spline');
% h111 = @(xa) h111r(xa) + 1i*h111i(xa);
% 
% h102r = griddedInterpolant(tcalc(:),(sol_1r(4,:)).','spline');
% h102i = griddedInterpolant(tcalc(:),(sol_1i(4,:)).','spline');
% h102 = @(xa) h102r(xa) + 1i*h102i(xa);
% 
% h012r = griddedInterpolant(tcalc(:),(sol_1r(5,:)).','spline');
% h012i = griddedInterpolant(tcalc(:),(sol_1i(5,:)).','spline');
% h012 = @(xa) h012r(xa) + 1i*h012i(xa);
% 
% h300 = (ap^3*gamma*L3)/(lf - 3*ls);
% h030 = (conj(ap)^3*gamma*L3)/(lf - 3*conj(ls)); 
% h210 = -(3* ap^2 *gamma*L3*conj(ap))/(-lf + 2*ls + conj(ls)); 
% h120 = -(3*ap*gamma*L3*conj(ap)^2)/(-lf + ls + 2*conj(ls));
% 
% 
% L_0 = -3000;
% tspan_0 = linspace(L_0,300,N_Max);
% positive_Lvec = find(tspan_0>0);
% positive_L = positive_Lvec(1)-1;
% 
% 
% h_coeff_1=@(t,ta) [exp((lf-3*ls)*(t-ta)).*(gamma.*L1.*h102(ta).*ap^3 +gamma.*conj(L1).*h012(ta).*ap.^3 - gamma.*L3.*3.*xi_sol1(ta).^2.*(bp.*h300+conj(bp.*h030)) - gamma.*L3.*3.*(ap.^2.*(bp.*h102(ta) + conj(bp.*h012(ta)) ) ) + h111(ta).*conj(L1).*gamma.*3.*xi_sol1(ta).*ap^2 + h201(ta).*6.*xi_sol1(ta).*ap.^2.*L1.*gamma );
% exp((lf-3*conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*conj(ap)^3 +gamma*conj(L1)*h012(ta)*conj(ap)^3 - gamma*L3*3*xi_sol1(ta).^2.*(bp*h030+conj(bp*h300)) - gamma*L3*3*(conj(ap)^2*(bp*h012(ta) + conj(bp*h102(ta)) ) ) + h111(ta).*(L1).*gamma.*3.*xi_sol1(ta).*conj(ap)^2 + h021(ta).*6.*xi_sol1(ta).*conj(ap)^2.*conj(L1).*gamma );
% exp((lf-ls-2*conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*3*ap*conj(ap)^2 + gamma*conj(L1)*h102(ta)*3*ap*conj(ap)^2-gamma*L3*3*xi_sol1(ta).^2.*(bp*h120+conj(bp*h210)) - gamma*L3*3*(conj(ap)^2*(bp*h102(ta)+conj(bp*h012(ta))) + 2*ap*conj(ap)*(bp*h012(ta)+conj(bp*h102(ta))) ) + h111(ta).*(L1)*gamma*3.*xi_sol1(ta)*ap*conj(ap)*2 + h021(ta).*6.*xi_sol1(ta)*2*ap*conj(ap)*conj(L1)*gamma );
% exp((lf-2*ls-conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*3*ap^2*conj(ap) + gamma*conj(L1)*h102(ta)*3*ap^2*conj(ap)-gamma*L3*3*xi_sol1(ta).^2.*(bp*h210+conj(bp*h120)) - gamma*L3*3*((ap)^2*(bp*h012(ta)+conj(bp*h102(ta))) + 2*ap*conj(ap)*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*conj(L1)*gamma*3.*xi_sol1(ta)*ap*conj(ap)*2 + h201(ta).*6.*xi_sol1(ta)*2*ap*conj(ap)*(L1)*gamma );
% exp((lf-4*ls)*(t-ta)).*(gamma*conj(L1)*ap^3*h111(ta) + 2*L1*gamma*ap^3*h201(ta) - gamma*L3*3*ap^2*(bp*h201(ta)+conj(bp*h021(ta))) - gamma*L3*3*xi_sol1(ta).*2.*(ap*(bp*h300 + conj(bp*h030))) );
% exp((lf-4*conj(ls))*(t-ta)).*(gamma*L1*conj(ap)^3*h111(ta) + 2* gamma *conj(L1)*conj(ap)^3*h021(ta) -gamma*L3*3*conj(ap)^2*(bp*h021(ta)+conj(bp*h201(ta))) - gamma*L3*3*xi_sol1(ta).*2*(conj(ap)*(bp*h030+conj(bp*h300))) );
% exp((lf-3*conj(ls)-ls)*(t-ta)).*(h111(ta).*(gamma*L1*3*conj(ap)^2*ap+gamma*conj(L1)*conj(ap)^3) + 2*h201(ta)*gamma*L1*conj(ap)^3 + 2*h021(ta)*gamma*conj(L1)*3*conj(ap)^2*ap -3*gamma*L3*(conj(ap)^2*(bp*h111(ta)+conj(bp*h111(ta))) + 2*ap*conj(ap)*(bp*h021(ta) + conj(bp*h201(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h030+ap*conj(bp*h300)+conj(ap)*bp*h120+conj(ap)*conj(bp*h210)) );
% exp((lf-conj(ls)-3*ls)*(t-ta)).*(h111(ta).*(gamma*conj(L1)*3*conj(ap)*ap^2+gamma*L1*(ap)^3) + 2*h021(ta)*gamma*conj(L1)*(ap)^3 + 2*h201(ta)*gamma*(L1)*3*(ap)^2*conj(ap) -3*gamma*L3*((ap)^2*(bp*h111(ta)+conj(bp*h111(ta))) + 2*ap*conj(ap)*(bp*h201(ta) + conj(bp*h021(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h210+ap*conj(bp*h120)+conj(ap)*bp*h300+conj(ap)*conj(bp*030))  );
% exp((lf-2*conj(ls)-2*ls)*(t-ta)).*(h111(ta).*(L1*gamma*3*ap^2*conj(ap)+conj(L1)*gamma*3*ap*conj(ap)^2) + 2*h021(ta)*gamma*conj(L1)*3*ap^2*conj(ap) + 2*h201(ta)*gamma*(L1)*3*(ap)*conj(ap)^2 -3*gamma*L3*((ap)^2*(bp*h021(ta)+conj(bp*h201(ta))) + (conj(ap))^2*(bp*h201(ta)+conj(bp*h021(ta))) + 2*ap*conj(ap)*(bp*h111(ta) + conj(bp*h111(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h120+ap*conj(bp*h210)+conj(ap)*bp*h210+conj(ap)*conj(bp*120))  );
% exp((lf-ls-conj(ls))*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*2*ap*conj(ap) -L3*gamma*3*xi_sol1(ta).^2.*(bp*h111(ta)+conj(bp*h111(ta))) - gamma*L3*6*xi_sol1(ta).*(ap*(bp*h012(ta)+conj(bp*h102(ta))) + conj(ap)*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*L1*gamma*3.*xi_sol1(ta).^2*conj(ap)+ conj(L1)*gamma*3*xi_sol1(ta).^2*ap.*h111(ta) + h201(ta).*2.*L1*gamma*3.*xi_sol1(ta).^2*conj(ap) + h021(ta).*2*L1*gamma*3*ap.*xi_sol1(ta).^2 + h102(ta).*gamma*L1*3*2*ap*conj(ap).*xi_sol1(ta) + h012(ta).*gamma*conj(L1)*3*2*ap*conj(ap).*xi_sol1(ta) );
% exp((lf-2*ls)*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*ap^2 -L3*gamma*3*xi_sol1(ta).^2.*(bp*h201(ta)+conj(bp*h021(ta))) - gamma*L3*6*xi_sol1(ta).*(ap*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*conj(L1)*gamma*3.*xi_sol1(ta).^2*(ap) + h201(ta).*2.*L1*gamma*3.*xi_sol1(ta).^2*(ap) + h102(ta).*gamma*L1*3.*xi_sol1(ta)*ap^2 + h012(ta).*gamma*conj(L1)*3.*xi_sol1(ta)*ap^2 );
% exp((lf-2*conj(ls))*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*conj(ap)^2 -L3*gamma*3*xi_sol1(ta).^2.*(bp*h021(ta)+conj(bp*h201(ta))) - gamma*L3*6*xi_sol1(ta).*(conj(ap)*(bp*h012(ta)+conj(bp*h102(ta))) ) + h111(ta).*(L1)*gamma*3.*xi_sol1(ta).^2*conj(ap) + h021(ta).*2.*conj(L1)*gamma*3.*xi_sol1(ta).^2*conj(ap) + h102(ta).*gamma*L1*3.*xi_sol1(ta)*conj(ap)^2 + h012(ta).*gamma*conj(L1)*3.*xi_sol1(ta)*conj(ap)^2 );
% exp((lf-ls)*(t-ta)).*(-gamma*L3*3*xi_sol1(ta).^2.*(bp*h102(ta)+conj(bp*h012(ta))) + gamma*L1*h102(ta).*3.*xi_sol1(ta).^2*ap + gamma*conj(L1)*h012(ta)*3.*xi_sol1(ta).^2*ap  );
% exp((lf-conj(ls))*(t-ta)).*(-gamma*L3*3*xi_sol1(ta).^2.*(bp*h012(ta)+conj(bp*h102(ta))) + gamma*L1*h102(ta).*3.*xi_sol1(ta).^2*conj(ap) + gamma*conj(L1)*h012(ta)*3.*xi_sol1(ta).^2*conj(ap))
%  ];

% exp((lf-3*ls)*(t-ta)).*(gamma*L1*h102(ta)*ap^3 +gamma*conj(L1)*h012(ta)*ap^3 - gamma*L3*3*xi_sol1(ta).^2.*(bp*h300+conj(bp*h030)) - gamma*L3*3*(ap^2*(bp*h102(ta) + conj(bp*h012(ta)) ) ) + h111(ta).*conj(L1)*gamma*3.*xi_sol1(ta)*ap^2 + h201(ta).*6.*xi_sol1(ta)*ap^2*L1*gamma )
% 
% exp((lf-3*conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*conj(ap)^3 +gamma*conj(L1)*h012(ta)*conj(ap)^3 - gamma*L3*3*xi_sol1(ta).^2.*(bp*h030+conj(bp*h300)) - gamma*L3*3*(conj(ap)^2*(bp*h012(ta) + conj(bp*h102(ta)) ) ) + h111(ta).*(L1)*gamma*3.*xi_sol1(ta)*conj(ap)^2 + h021(ta).*6.*xi_sol1(ta)*conj(ap)^2*conj(L1)*gamma )
% 
% exp((lf-ls-2*conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*3*ap*conj(ap)^2 + gamma*conj(L1)*h102(ta)*3*ap*conj(ap)^2-gamma*L3*3*xi_sol1(ta).^2.*(bp*h120+conj(bp*h210)) - gamma*L3*3*(conj(ap)^2*(bp*h102(ta)+conj(bp*h012(ta))) + 2*ap*conj(ap)*(bp*h012(ta)+conj(bp*h102(ta))) ) + h111(ta).*(L1)*gamma*3.*xi_sol1(ta)*ap*conj(ap)*2 + h021(ta).*6.*xi_sol1(ta)*2*ap*conj(ap)*conj(L1)*gamma )
% 
% exp((lf-2*ls-conj(ls))*(t-ta)).*(gamma*L1*h102(ta)*3*ap^2*conj(ap) + gamma*conj(L1)*h102(ta)*3*ap^2*conj(ap)-gamma*L3*3*xi_sol1(ta).^2.*(bp*h210+conj(bp*h120)) - gamma*L3*3*((ap)^2*(bp*h012(ta)+conj(bp*h102(ta))) + 2*ap*conj(ap)*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*conj(L1)*gamma*3.*xi_sol1(ta)*ap*conj(ap)*2 + h201(ta).*6.*xi_sol1(ta)*2*ap*conj(ap)*(L1)*gamma )

% % %h401 exp((lf-4*ls)*(t-ta)).*(gamma*conj(L1)*ap^3*h111(ta) + 2*L1*gamma*ap^3*h201(ta) - gamma*L3*3*ap^2*(bp*h201(ta)+conj(bp*h021(ta))) - gamma*L3*3*xi_sol1(ta).*2.*(ap*(bp*h300 + conj(bp*h030))) )
% % % 
% % %h041 exp((lf-4*conj(ls))*(t-ta)).*(gamma*L1*conj(ap)^3*h111(ta) + 2* gamma *conj(L1)*conj(ap)^3*h021(ta) -gamma*L3*3*conj(ap)^2*(bp*h021(ta)+conj(bp*h201(ta))) - gamma*L3*3*xi_sol1(ta).*2*(conj(ap)*(bp*h030+conj(bp*h300))) )
% % % 
% % %h131 exp((lf-3*conj(ls)-ls)*(t-ta)).*(h111(ta).*(gamma*L1*3*conj(ap)^2*ap+gamma*conj(L1)*conj(ap)^3) + 2*h201(ta)*gamma*L1*conj(ap)^3 + 2*h021(ta)*gamma*conj(L1)*3*conj(ap)^2*ap -3*gamma*L3*(conj(ap)^2*(bp*h111(ta)+conj(bp*h111(ta))) + 2*ap*conj(ap)*(bp*h021(ta) + conj(bp*h201(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h030+ap*conj(bp*h300)+conj(ap)*bp*h120+conj(ap)*conj(bp*h210)) )
% % % 
% % %h311 exp((lf-conj(ls)-3*ls)*(t-ta)).*(h111(ta).*(gamma*conj(L1)*3*conj(ap)*ap^2+gamma*L1*(ap)^3) + 2*h021(ta)*gamma*conj(L1)*(ap)^3 + 2*h201(ta)*gamma*(L1)*3*(ap)^2*conj(ap) -3*gamma*L3*((ap)^2*(bp*h111(ta)+conj(bp*h111(ta))) + 2*ap*conj(ap)*(bp*h201(ta) + conj(bp*h021(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h210+ap*conj(bp*h120)+conj(ap)*bp*h300+conj(ap)*conj(bp*030))  )
% % % 
% % %h221 exp((lf-2*conj(ls)-2*ls)*(t-ta)).*(h111(ta).*(L1*gamma*3*ap^2*conj(ap)+conj(L1)*gamma*3*ap*conj(ap)^2) + 2*h021(ta)*gamma*conj(L1)*3*ap^2*conj(ap) + 2*h201(ta)*gamma*(L1)*3*(ap)*conj(ap)^2 -3*gamma*L3*((ap)^2*(bp*h021(ta)+conj(bp*h201(ta))) + (conj(ap))^2*(bp*h201(ta)+conj(bp*h021(ta))) + 2*ap*conj(ap)*(bp*h111(ta) + conj(bp*h111(ta))) ) - gamma*L3*3*2*xi_sol1(ta).*(ap*bp*h120+ap*conj(bp*h210)+conj(ap)*bp*h210+conj(ap)*conj(bp*120))  ) 

% %h113 exp((lf-ls-conj(ls))*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*2*ap*conj(ap) -L3*gamma*3*xi_sol1(ta).^2.*(bp*h111(ta)+conj(bp*h111(ta))) - gamma*L3*6*xi_sol1(ta).*(ap*(bp*h012(ta)+conj(bp*h102(ta))) + conj(ap)*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*L1*gamma*3*xi_sol1(ta).^2*conj(ap)+ conj(L1)*gamma*3*xi_sol1(ta).^2*ap.*h111(ta) + h201(ta).*2.*L1*gamma*3*xi_sol1(ta).^2*conj(ap) + h021(ta).*2*L1*gamma*3*ap*xi_sol1(ta).^2 + h102(ta).*gamma*L1*3*2*ap*conj(ap).*xi_sol1(ta) + h012(ta).*gamma*conj(L1)*3*2*ap*conj(ap).*xi_sol1(ta) )
% % 
% %h203 exp((lf-2*ls)*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*ap^2 -L3*gamma*3*xi_sol1(ta).^2.*(bp*h201(ta)+conj(bp*h021(ta))) - gamma*L3*6*xi_sol1(ta).*(ap*(bp*h102(ta)+conj(bp*h012(ta))) ) + h111(ta).*conj(L1)*gamma*3*xi_sol1(ta).^2*(ap) + h201(ta).*2.*L1*gamma*3*xi_sol1(ta).^2*(ap) + h102(ta).*gamma*L1*3.*xi_sol1(ta)*ap^2 + h012(ta).*gamma*conj(L1)*3.*xi_sol1(ta)*ap^2 )
% % 
% %h023 exp((lf-2*conj(ls))*(t-ta)).*(-gamma*L3*3*xi3_sol1(ta)*conj(ap)^2 -L3*gamma*3*xi_sol1(ta).^2.*(bp*h021(ta)+conj(bp*h201(ta))) - gamma*L3*6*xi_sol1(ta).*(conj(ap)*(bp*h012(ta)+conj(bp*h102(ta))) ) + h111(ta).*(L1)*gamma*3*xi_sol1(ta).^2*conj(ap) + h021(ta).*2.*conj(L1)*gamma*3*xi_sol1(ta).^2*conj(ap) + h102(ta).*gamma*L1*3.*xi_sol1(ta)*conj(ap)^2 + h012(ta).*gamma*conj(L1)*3.*xi_sol1(ta)*conj(ap)^2 )

% h104 exp((lf-ls)*(t-ta)).*(-gamma*L3*3*xi_sol1(ta).^2.*(bp*h102(ta)+conj(bp*h012(ta))) + gamma*L1*h102(ta).*3.*xi_sol1(ta).^2*ap + gamma*conj(L1)*h012(ta)*3.*xi_sol1(ta).^2*ap  )
% h014 exp((lf-conj(ls))*(t-ta)).*(-gamma*L3*3*xi_sol1(ta).^2.*(bp*h012(ta)+conj(bp*h102(ta))) + gamma*L1*h102(ta).*3.*xi_sol1(ta).^2*conj(ap) + gamma*conj(L1)*h012(ta)*3.*xi_sol1(ta).^2*conj(ap))


% % tcalc = tspan_0(positive_L:end);
% % 
% % sol_1i = [];
% % sol_1r = [];
% % tic
% % for i = 1:max(size(tcalc))
% % Lo = positive_L - 1 + i ;
% % T = tspan_0(start_time(1):Lo);
% % ti = tcalc(i);
% % 
% % INx = h_coeff_1(ti,T);
% % INxr = real(INx);
% % INxi = imag(INx);
% % Q1xr = trapz(T,INxr,2);
% % Q1xi = trapz(T,INxi,2);
% % sol_1r  =[sol_1r,Q1xr];
% % sol_1i  =[sol_1i,Q1xi];
% % end    
% % toc
% % 
% % h302r = griddedInterpolant(tcalc(:),(sol_1r(1,:)).','spline');
% % h302i = griddedInterpolant(tcalc(:),(sol_1i(1,:)).','spline');
% % h302 = @(xa) h302r(xa) + 1i*h302i(xa);
% % 
% % h032r = griddedInterpolant(tcalc(:),(sol_1r(2,:)).','spline');
% % h032i = griddedInterpolant(tcalc(:),(sol_1i(2,:)).','spline');
% % h032 = @(xa) h032r(xa) + 1i*h032i(xa);
% % 
% % h122r = griddedInterpolant(tcalc(:),(sol_1r(3,:)).','spline');
% % h122i = griddedInterpolant(tcalc(:),(sol_1i(3,:)).','spline');
% % h122 = @(xa) h122r(xa) + 1i*h122i(xa);
% % 
% % h212r = griddedInterpolant(tcalc(:),(sol_1r(4,:)).','spline');
% % h212i = griddedInterpolant(tcalc(:),(sol_1i(4,:)).','spline');
% % h212 = @(xa) h212r(xa) + 1i*h212i(xa);
% % 
% % h401r = griddedInterpolant(tcalc(:),(sol_1r(5,:)).','spline');
% % h401i = griddedInterpolant(tcalc(:),(sol_1i(5,:)).','spline');
% % h401 = @(xa) h401r(xa) + 1i*h401i(xa);
% % 
% % h041r = griddedInterpolant(tcalc(:),(sol_1r(6,:)).','spline');
% % h041i = griddedInterpolant(tcalc(:),(sol_1i(6,:)).','spline');
% % h041 = @(xa) h041r(xa) + 1i*h041i(xa);
% % 
% % h131r = griddedInterpolant(tcalc(:),(sol_1r(7,:)).','spline');
% % h131i = griddedInterpolant(tcalc(:),(sol_1i(7,:)).','spline');
% % h131 = @(xa) h131r(xa) + 1i*h131i(xa);
% % 
% % h311r = griddedInterpolant(tcalc(:),(sol_1r(8,:)).','spline');
% % h311i = griddedInterpolant(tcalc(:),(sol_1i(8,:)).','spline');
% % h311 = @(xa) h311r(xa) + 1i*h311i(xa);
% % 
% % h221r = griddedInterpolant(tcalc(:),(sol_1r(9,:)).','spline');
% % h221i = griddedInterpolant(tcalc(:),(sol_1i(9,:)).','spline');
% % h221 = @(xa) h221r(xa) + 1i*h221i(xa);
% % 
% % h113r = griddedInterpolant(tcalc(:),(sol_1r(10,:)).','spline');
% % h113i = griddedInterpolant(tcalc(:),(sol_1i(10,:)).','spline');
% % h113 = @(xa) h113r(xa) + 1i*h113i(xa);
% % 
% % h203r = griddedInterpolant(tcalc(:),(sol_1r(11,:)).','spline');
% % h203i = griddedInterpolant(tcalc(:),(sol_1i(11,:)).','spline');
% % h203 = @(xa) h203r(xa) + 1i*h203i(xa);
% % 
% % h023r = griddedInterpolant(tcalc(:),(sol_1r(12,:)).','spline');
% % h023i = griddedInterpolant(tcalc(:),(sol_1i(12,:)).','spline');
% % h023 = @(xa) h023r(xa) + 1i*h023i(xa);
% % 
% % h104r = griddedInterpolant(tcalc(:),(sol_1r(13,:)).','spline');
% % h104i = griddedInterpolant(tcalc(:),(sol_1i(13,:)).','spline');
% % h104 = @(xa) h104r(xa) + 1i*h104i(xa);
% % 
% % h014r = griddedInterpolant(tcalc(:),(sol_1r(14,:)).','spline');
% % h014i = griddedInterpolant(tcalc(:),(sol_1i(14,:)).','spline');
% % h014 = @(xa) h014r(xa) + 1i*h014i(xa);
% % 
% % save('Coeff_SSM_NF_Exact.mat','h014','h104','h023','h203','h113','h221','h311','h131','h041','h401','h212','h122','h032','h302','h012','h102','h201','h021','h111');
% % 

% h302 = griddedInterpolant(tcalc(:),(sol_1(1,:)).','spline');
% h032 = griddedInterpolant(tcalc(:),(sol_1(2,:)).','spline');
% h122 = griddedInterpolant(tcalc(:),(sol_1(3,:)).','spline');
% h212 = griddedInterpolant(tcalc(:),(sol_1(4,:)).','spline');


%%


% h300 = (ap^3*gamma*L3)/(lf - 3*ls);
% h030 = (conj(ap)^3*gamma*L3)/(lf - 3*conj(ls)); 
% h210 = -(3* ap^2 *gamma*L3*conj(ap))/(-lf + 2*ls + conj(ls)); 
% h120 = -(3*ap*gamma*L3*conj(ap)^2)/(-lf + ls + 2*conj(ls));
% load('Coeff_SSM_NF_Exact.mat')
% 
% Order = 5;
% 
% ROM_ODE_ns =@(t,zk) rom_temp_model_noscale(t,zk,ls,lf,V,inv(V),gamma,epsilon,Order,XI_hyperbolic_solution,h300,h030,h210,h120,h012,h021,h032,h102,h111,h122,h201,h212,h302,h014,h104,h023,h203,h113,h221,h311,h131,h041,h401);
% ctspan = linspace(0,200,5000);
% 
% q0 = 0.8*exp(1i*0.3);
% y0 =compute_SSM_phy(XI_hyperbolic_solution,0,lf,ls,gamma,V,inv(V),q0,h300,h030,h210,h120,h012,h021,h032,h102,h111,h122,h201,h212,h302,h014,h104,h023,h203,h113,h221,h311,h131,h041,h401,epsilon);
% 
% IC = [real(q0);imag(q0)];
% 
% tic 
% [t_sol,ROM_sol] = ode45(ROM_ODE_ns,ctspan,IC);
% toc
% 
% 
% [y,modal] = compute_SSM_phy(XI_hyperbolic_solution,t_sol.',lf,ls,gamma,V,inv(V),(ROM_sol(:,1)+1i*ROM_sol(:,2)).',h300,h030,h210,h120,h012,h021,h032,h102,h111,h122,h201,h212,h302,h014,h104,h023,h203,h113,h221,h311,h131,h041,h401,epsilon);
% 
% tic
%  [F, lambda, V1, G, DG] = functionFromTensors(M, C, K, fnl);
% toc
% 
% off = 0;
% Gh = @(t,x) F(t,x)+epsilon*[0;0;Force_Lorenz(t);Force_Lorenz(t)];
% Solution = epsilon*XI_hyperbolic_solution{1,1}(t_sol.');
% 
% tic
% IC = y0;
% [tSP,SP_Traj] = ode45(Gh,ctspan,IC);
% toc
% 
% 
% figure
% indexR = 2;
% % plot(tSP,SP_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0,1])
% % hold on 
% % plot(t_sol,y(indexR,:),'--','LineWidth',3,'color',[1, 0, 0, 1])
% % hold on 
% plot(t_sol,Solution(:,indexR),'-','LineWidth',3,'color',[0,1,0,0.3])
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$q_2$','Interpreter','latex');
% % legend('FOM','ROM - O(\epsilon^{k_1} u^{k_2}, k_1 + k_2 =5)','O(\epsilon^5): x_{\epsilon}(t) Anchor Trajectory')
% legend('O(\epsilon^5): x_{\epsilon}(t) Anchor Trajectory')
% title('\epsilon = 0.1 and |\epsilon F_{max}| = 10')
% axis([0 100 -0.2 0.2])
% 
% tc = 0:1/1e3:2;
% yc = chirp(t,0,1,250);


% axes('position',[.25 .175 .65 .25])
% box on % put box around new pair of axes
% indexR = 1;
% stick = linspace(max(abs(SP_Traj(:,indexR))),-max(abs(SP_Traj(:,indexR))),200);
% 
% plot(t1f*ones(1,max(size(stick))),stick,'-','LineWidth',3,'color',[0.4940 0.1840 0.5560])
% hold on
% plot(t2f*ones(1,max(size(stick))),stick,'-','LineWidth',3,'color',[0.4940 0.1840 0.5560])
% hold on 
% plot(t3f*ones(1,max(size(stick))),stick,'-','LineWidth',3,'color',[0.4940 0.1840 0.5560])
% 
% hold on
% 
% plot(t_sol(k1enter),Solution(k1enter,indexR),'-o','LineWidth',1,'color','green')
% hold on
% plot(t_sol(k2enter),Solution(k2enter,indexR),'-o','LineWidth',1,'color','green')
% hold on 
% plot(t_sol(k3enter),Solution(k3enter,indexR),'-o','LineWidth',1,'color','green')
% hold on 
% plot(t_sol(k4enter),Solution(k4enter,indexR),'-o','LineWidth',1,'color','green')
% hold on 
% 
% plot(tSP(k1enter),SP_Traj(k1enter,indexR),'-','LineWidth',3,'color','black')
% hold on
% plot(t_sol(k1enter),y(indexR,k1enter),'--','LineWidth',3,'color','red')
% hold on 
% 
% plot(tSP(k2enter),SP_Traj(k2enter,indexR),'-','LineWidth',3,'color','black')
% hold on 
% plot(tSP(k3enter),SP_Traj(k3enter,indexR),'-','LineWidth',3,'color','black')
% hold on 
% plot(tSP(k4enter),SP_Traj(k4enter,indexR),'-','LineWidth',3,'color','black')
% hold on 
% 
% 
% 
%  
% plot(t_sol(k2enter),y(indexR,k2enter),'--','LineWidth',3,'color','red')
% hold on 
% plot(t_sol(k3enter),y(indexR,k3enter),'--','LineWidth',3,'color','red')
% hold on 
% plot(t_sol(k4enter),y(indexR,k4enter),'--','LineWidth',3,'color','red')
% hold on 
% 
% axis([41 160 -0.06 0.06])



% figure
% H = h111(ctspan);
% subplot(2,1,1)
% plot(tSP,real(H),'-','LineWidth',3,'color','black')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$Re[h^{111}(t)]$','Interpreter','latex');
% title('h^{111}(t)')
% subplot(2,1,2)
% plot(tSP,imag(H),'-','LineWidth',3,'color','black')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$Im[h^{111}(t)]$','Interpreter','latex');


% subplot(2,2,2)
% plot(tSP,SP_Traj(:,2),'-','LineWidth',3,'color','black')
% hold on 
% plot(tSP_ns,SP_Traj_ns(:,2),'-','LineWidth',3,'color','magenta')
% hold on
% plot(t_sol,y(2,:),'--','LineWidth',3,'color','red')
% hold on 
% plot(t_sol,y_ns(2,:),'--','LineWidth',3,'color','blue')
% hold on
% plot(t_sol,Solution(2,:),'-','LineWidth',3,'color','green')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$q_2 \,[$m$]$','Interpreter','latex');
% legend('FOM','FOM (no scaling)','ROM','ROM (no scaling)','$O(\epsilon)$ Solution')
% 
% subplot(2,2,3)
% plot(tSP,SP_Traj(:,3),'-','LineWidth',3,'color','black')
% hold on 
% plot(tSP_ns,SP_Traj_ns(:,3),'-','LineWidth',3,'color','magenta')
% hold on
% plot(t_sol,y(3,:),'--','LineWidth',3,'color','red')
% hold on 
% plot(t_sol,y_ns(3,:),'--','LineWidth',3,'color','blue')
% hold on
% plot(t_sol,Solution(3,:),'-','LineWidth',3,'color','green')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$\dot{q}_1 \,[$m/s$]$','Interpreter','latex');
% legend('FOM','FOM (no scaling)','ROM','ROM (no scaling)','$O(\epsilon)$ Solution')
% 
% subplot(2,2,4)
% plot(tSP,SP_Traj(:,4),'-','LineWidth',3,'color','black')
% hold on 
% plot(tSP_ns,SP_Traj_ns(:,4),'-','LineWidth',3,'color','magenta')
% hold on
% plot(t_sol,y(4,:),'--','LineWidth',3,'color','red')
% hold on 
% plot(t_sol,y_ns(4,:),'--','LineWidth',3,'color','blue')
% hold on
% plot(t_sol,Solution(4,:),'-','LineWidth',3,'color','green')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$\dot{q}_2 \,[$m/s$]$','Interpreter','latex');
% legend('FOM','FOM (no scaling)','ROM','ROM (no scaling)','$O(\epsilon)$ Solution')
% 
% 
% sgtitle('\epsilon = 1')
% 
% figure 
% plot3(SP_Traj(:,1),SP_Traj(:,2),SP_Traj(:,3),'--','LineWidth',3,'color','black')
% hold on 
% plot3(Solution(1,:),Solution(2,:),Solution(3,:),'-','LineWidth',3,'color','red')
% 
% 

% Figure
%%
% hFig = figure('DefaultAxesFontSize',20);                       % Bring up new figure
% imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
% hFig.WindowState = 'maximized';
% pause(0.5)
% delta = 0.5;
% tl = 800;
% clear movieVector
% shift_traj = y;
% shift_trajZ = [real(modal(1,:));imag(modal(1,:));imag(modal(4,:));real(modal(4,:))];%y;
% ind = 1;
% SP_TrajZ = inv(V)*SP_Traj.';
% SP_TrajZ = [real(SP_TrajZ(1,:));imag(SP_TrajZ(1,:));imag(SP_TrajZ(4,:));real(SP_TrajZ(4,:))];
% SolutionZ = inv(V)*Solution;
% SolutionZ = [real(SolutionZ(1,:));imag(SolutionZ(1,:));imag(SolutionZ(4,:));real(SolutionZ(4,:))];
% for i = 4200:4200%(max(size(t_sol))-tl)
% 
%      clf
% subplot(2, 2, [1 3])
%      rhosamp = linspace(0,2,201);
%      plotdofs = [1 2 4]; 
%      [RHO,Theta]=meshgrid(linspace(0,2,201),linspace(0,2*pi,201));
%      Z = epsilon^(1/Order)*RHO.*exp(1i*Theta);
%      t_ssm = ones(1,201*201)*t_sol(i+tl);
%      [ZZ,modalZ] = compute_SSM_phy(xi_net_sol,t_ssm,lf,ls,gamma,V,inv(V),Z(:).');
%      %Z1 = reshape(real(modalZ(1,:)),201,201);
%      %Z2 = reshape(imag(modalZ(1,:)),201,201);
%      %Z3 = reshape(real(modalZ(4,:)),201,201);
%      Z1 = reshape(ZZ(1,:),201,201);
%      Z2 = reshape(ZZ(2,:),201,201);
%      Z3 = reshape(ZZ(3,:),201,201);
% 
%      h = surf(Z1,Z2,Z3)
%      h.EdgeColor = 'none';
%      h.FaceColor = 'green';
%      h.FaceAlpha = 0.3;
%      hold on 
% 
% %      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
% %      hold on 
%      starto = 4200;
%      plot3(Solution(plotdofs(1),starto:end),Solution(plotdofs(2),starto:end),Solution(plotdofs(3),starto:end),'--','LineWidth',3,'color','red')
%      hold on 
%      plot3(Solution(plotdofs(1),starto),Solution(plotdofs(2),starto),Solution(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(Solution(plotdofs(1),end),Solution(plotdofs(2),end),Solution(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
% 
%      plot3(shift_traj(plotdofs(1),i:i+tl),shift_traj(plotdofs(2),i:i+tl),shift_traj(plotdofs(3),i:i+tl),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
%      hold on
%      plot3(shift_traj(plotdofs(1),i+tl),shift_traj(plotdofs(2),i+tl),shift_traj(plotdofs(3),i+tl),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
%      hold on
% 
% %      plot3(shift_trajZ(plotdofs(1),i:i+tl),shift_trajZ(plotdofs(2),i:i+tl),shift_trajZ(plotdofs(3),i:i+tl),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %      hold on
% %      plot3(shift_trajZ(plotdofs(1),i+tl),shift_trajZ(plotdofs(2),i+tl),shift_trajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %      hold on 
% %      plot3(SP_TrajZ(plotdofs(1),i:i+tl),SP_TrajZ(plotdofs(2),i:i+tl),SP_TrajZ(plotdofs(3),i:i+tl),'-','LineWidth',3,'Color','black');
% %      hold on
% %      plot3(SP_TrajZ(plotdofs(1),i+tl),SP_TrajZ(plotdofs(2),i+tl),SP_TrajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color','black');
% %     
% %      plot3(SP_Traj(i:i+tl,plotdofs(1)),SP_Traj(i:i+tl,plotdofs(2)),SP_Traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color','black');
% %      hold on
% %      plot3(SP_Traj(i+tl,plotdofs(1)),SP_Traj(i+tl,plotdofs(2)),SP_Traj(i+tl,plotdofs(3)),'.','MarkerSize',30,'Color','black');
% %       hold on 
%      
%     legend('SSM - $W_{\epsilon}(x_{\epsilon}(t))$','$O(\epsilon)$ solution','Start','End','ROM','ROM start')%,'FOM','FOM head')
%  
% %     h
% %     hold on
% 
% 
%      xlabel('$q_1 [m]$','Interpreter','latex');
%      ylabel('$q_2 [m]$','Interpreter','latex');
%      zlabel('$\dot{q}_2 [m/s]$','Interpreter','latex');
%         daspect([1,1,1])
% 
%     view(32,26)
%     grid on
% %      W0 = WTF(target(Omega*traj.time(i+tl)));
% %   
% % coord_change_vec2 = 0;
% % coord_change_vec = [traj_actual(i+tl,plotdofs(1));traj_actual(i+tl,plotdofs(2));traj_actual(i+tl,plotdofs(3))];
% %      [zdof1,zdof2,zdof3]=plot_2D_auto_SSM(coord_change_vec2,coord_change_vec,W0,rhosamp,plotdofs,{'$q_1$ m','$q_2$ m','$\dot{q}_2$ m/s'});
% %      hold on 
% %      plot3(traj_actual1(:,plotdofs(1)),traj_actual1(:,plotdofs(2)),traj_actual1(:,plotdofs(3)),'-','Linewidth',3,'color','red');
% %     hold on
% %     plot3(coord_change_vec(1),coord_change_vec(2),coord_change_vec(3),'.','MarkerSize',30,'color','blue');
% %     hold on
% %     plot3(shift_traj(i:i+tl,plotdofs(1)),shift_traj(i:i+tl,plotdofs(2)),shift_traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(shift_traj(i+tl,plotdofs(1)),shift_traj(i+tl,plotdofs(2)),shift_traj(i+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(shift_traj(i:i+tl,plotdofs(1)),shift_traj(i:i+tl,plotdofs(2)),shift_traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(SP_Traj_Off(i:i+2*tl,plotdofs(1)),SP_Traj_Off(i:i+2*tl,plotdofs(2)),SP_Traj_Off(i:i+2*tl,plotdofs(3)),'-','LineWidth',3,'Color','black');
% %     hold on 
% %     plot3(SP_Traj_Off(i+2*tl,plotdofs(1)),SP_Traj_Off(i+2*tl,plotdofs(2)),SP_Traj_Off(i+2*tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %     
% %     hold on 
% %     view([108 8])
% %     axis([-0.02,0.02,-0.02,0.02,-0.02,0.02])
% %    daspect([1 1 1])
%     axis([-0.1,0.1,-0.1,0.1,-0.4,0.4])
%     
% 
% subplot(2, 2, 2)
% rhosamp = linspace(0,2,201);
%      plotdofs = [1 2 4]; 
%      [RHO,Theta]=meshgrid(linspace(0,2,201),linspace(0,2*pi,201));
%      Z = epsilon^(1/Order)*RHO.*exp(1i*Theta);
%      t_ssm = ones(1,201*201)*t_sol(i+tl);
%      [ZZ,modalZ] = compute_SSM_phy(xi_net_sol,t_ssm,lf,ls,gamma,V,inv(V),Z(:).');
%      %Z1 = reshape(real(modalZ(1,:)),201,201);
%      %Z2 = reshape(imag(modalZ(1,:)),201,201);
%      %Z3 = reshape(real(modalZ(4,:)),201,201);
%      Z1 = reshape(ZZ(1,:),201,201);
%      Z2 = reshape(ZZ(2,:),201,201);
%      Z3 = reshape(ZZ(3,:),201,201);
% 
%      h = surf(Z1,Z2,Z3)
%      h.EdgeColor = 'none';
%      h.FaceColor = 'green';
%      h.FaceAlpha = 0.3;
%      hold on 
% 
% %      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
% %      hold on 
%      starto = 4000;
%      plot3(Solution(plotdofs(1),starto:end),Solution(plotdofs(2),starto:end),Solution(plotdofs(3),starto:end),'-','LineWidth',3,'color','red')
%      hold on 
%      plot3(Solution(plotdofs(1),starto),Solution(plotdofs(2),starto),Solution(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(Solution(plotdofs(1),end),Solution(plotdofs(2),end),Solution(plotdofs(3),end),'.','MarkerSize',30,'color',[0.4660 0.6740 0.1880])
%      hold on 
%      plot3(shift_traj(plotdofs(1),1:300),shift_traj(plotdofs(2),1:300),shift_traj(plotdofs(3),1:300),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
%      hold on
%      plot3(shift_traj(plotdofs(1),300),shift_traj(plotdofs(2),300),shift_traj(plotdofs(3),300),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
%      hold on
%      xlabel('$q_1 [m]$','Interpreter','latex');
%      ylabel('$q_2 [m]$','Interpreter','latex');
%      zlabel('$\dot{q}_2 [m/s]$','Interpreter','latex');
%         daspect([1,1,1])
% 
% subplot(2, 2, 4)
% rhosamp = linspace(0,2,201);
%      plotdofs = [1 2 4]; 
%      [RHO,Theta]=meshgrid(linspace(0,2,201),linspace(0,2*pi,201));
%      Z = epsilon^(1/Order)*RHO.*exp(1i*Theta);
%      t_ssm = ones(1,201*201)*t_sol(i+tl);
%      [ZZ,modalZ] = compute_SSM_phy(xi_net_sol,t_ssm,lf,ls,gamma,V,inv(V),Z(:).');
%      Z1 = reshape(real(modalZ(1,:)),201,201);
%      Z2 = reshape(imag(modalZ(1,:)),201,201);
%      Z3 = reshape(real(modalZ(4,:)),201,201);
% %      Z1 = reshape(ZZ(1,:),201,201);
% %      Z2 = reshape(ZZ(2,:),201,201);
% %      Z3 = reshape(ZZ(3,:),201,201);
% 
%      h = surf(Z1,Z2,Z3)
%      h.EdgeColor = 'none';
%      h.FaceColor = 'green';
%      h.FaceAlpha = 0.3;
%      hold on 
% 
% %      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
% %      hold on 
%      starto = 4000;
%      plot3(SolutionZ(plotdofs(1),starto:end),SolutionZ(plotdofs(2),starto:end),SolutionZ(plotdofs(3),starto:end),'-','LineWidth',3,'color','red')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),starto),SolutionZ(plotdofs(2),starto),SolutionZ(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color',[0.4660 0.6740 0.1880])
%     hold on
%      plot3(shift_trajZ(plotdofs(1),1:300),shift_trajZ(plotdofs(2),1:300),shift_trajZ(plotdofs(3),1:300),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
%      hold on
%      plot3(shift_trajZ(plotdofs(1),300),shift_trajZ(plotdofs(2),300),shift_trajZ(plotdofs(3),300),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
%      hold on
%     xlabel('$Re(z)$','Interpreter','latex');
%      ylabel('$Imag(z)$','Interpreter','latex');
%      zlabel('$Re(p)$','Interpreter','latex');
%         daspect([1,1,1])
% 
% set(gcf,'color','white')
%     figssm = gcf;
% 
% movieVector(ind) = getframe(hFig);
%     
%     ind = ind +1;
% 
% end

%%


% figure 
% plot(tSP,SP_TrajZ(4,:),'-','LineWidth',3,'color','black')
% hold on 
% %plot(tSP_ns,SP_Traj_ns(:,1),'-','LineWidth',3,'color','magenta')
%   plot(t_sol,shift_trajZ(4,:),'--','LineWidth',3,'color','red')
%   hold on 
% %plot(t_sol,y_ns(1,:),'--','LineWidth',3,'color','blue')
%   hold on
%  plot(t_sol,SolutionZ(4,:),'--','LineWidth',3,'color','green')
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$q_1 \,[$m$]$','Interpreter','latex');
% legend('FOM','FOM (no scaling)','ROM','ROM (no scaling)','$O(\epsilon)$ Solution')

% myWriter = VideoWriter('SSM_Shaw_Pierre_Chaotic_shaking_slow','MPEG-4');
% myWriter.FrameRate = 10;
% 
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);

%Processing figures in Matlab
% origUnits = fig.Units;
% fig.Units = fig.PaperUnits;
% fig.PaperSize = fig.Position(3:4);
% fig.Units = origUnits;
% exportgraphics(fig, 'output.pdf');


% 
% % DF_a_new = 2*epsilon.*(sigma*(yy(:,1).*(rho - yy(:,3))-yy(:,2)-sigma.*(yy(:,2)-yy(:,1)))/(max(yy(:,1))));
% 
% % [dXdFPa, XtruncPa,F_truncPa] = ftd(F_a_new.', tspan);
% %  
% % DForce_Lorenz = @(xq) interp1(F_truncPa(:),dXdFPa(:),xq,'spline');
% 
% 
% 
% % Load_full_a = zeros(n,1000);
% % Load_full_a(2:3:n,:) = repmat(F_a_new.',nElements,1);
% % 
% % figure
% % plot(tk,F_a_new)
% % 
% % figure 
% % plot3(yy(:,1),yy(:,2),yy(:,3));
% 
% 
% %% Quasi Static Response for transverse loads
% 
% 
% % disp('Starting Quasi-Static response estimator') % 33 minutes roughly to compute the estimator
% % % % tic
% % % IC = getStaticResponse(K, M, F, Load_full_a, 1, PlotFieldonDefMesh);
% % % toc
% % % ICL = IC(1:n,:);
% % % save('Quasi_Static_Lorenz.mat','ICL')
% % % load('Quasi_Static.mat','IC')
% % 
% 
% % N_grid = 1000;
% % Grid_F = linspace(-max(F_a_new),max(F_a_new),N_grid);
% % WT = cell(1,N_grid);
% % RT = cell(1,N_grid);
% % DJ = cell(1,N_grid);
% % ICTT = cell(1,N_grid);
% % spec_ratio = [];
% % backbonecurve = cell(1,N_grid);
% % for i = 1:N_grid
% % load_vector = zeros(n,1);
% % load_vector(2:3:n,1) = Grid_F(i);
% % IC = getStaticResponse(K, M, F, load_vector, 0, 0);
% % ICTT{1,i} = IC(1:n,1);
% % [K_shift,fnl_shift] = build_model_shift(K,fnl,IC(1:n,1));
% % DS = DynamicalSystem();
% % set(DS,'M',M,'C',C,'K',K_shift,'fnl',fnl_shift);
% % set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% % [V,D,W] = DS.linear_spectral_analysis();
% % V = reorient_vec(V);
% % W = reorient_vec(W);
% % DS.spectrum.V = V;
% % DS.spectrum.W = W;
% % 
% % spec_ratio = [spec_ratio,real(D(3))/real(D(1))];
% % S = SSM(DS);
% % set(S.Options, 'reltol', 1,'notation','multiindex');
%   resonant_modes = [1 2]; % choose master spectral subspace
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
% % DJ{1,i} = lamdMaster;
% % 
% % % hj = matlabFunction(y(2));
% % % backbonecurve{1,i} = hj;
% % end    
% % 
% %  save('RTSP.mat','RT','R0','Grid_F');
% %  save('WTSP.mat','WT','W0','Grid_F');
% %  save('DSP.mat','DJ','Grid_F');
% %  save('ICSP.mat','ICTT','Grid_F');
% 
% 
% %  lambda = [];
% %  for i = 1:1000
% %     lambda = [lambda,real(DJ{1,i}(1,1))];
% %  end    
% 
% 
% % X_s = IC(1:n,1);
% % [K_shift,fnl_shift] = build_model_shift(K,fnl,X_s);
% % DS = DynamicalSystem();
% % set(DS,'M',M,'C',C,'K',K_shift,'fnl',fnl_shift);
% % set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% % [V,D1,W] = DS.linear_spectral_analysis();
% 
% %%
% 
% % load('RTbeam4.mat','RT','R0','Grid_F');
% % load('WTbeam4.mat','WT','W0','Grid_F');
% % load('Dbeam4.mat','DJ','Grid_F');
% % load('ICbeam4.mat','ICTT','Grid_F');
% 
% % load('RTSP.mat','RT','R0','Grid_F');
% % load('WTSP.mat','WT','W0','Grid_F');
% % load('DSP.mat','DJ','Grid_F');
% % load('ICSP.mat','ICTT','Grid_F');
% % 
% % RF = grid_SSM_interp1(RT,Grid_F,N_grid);
% % WF = grid_SSM_interp1(WT,Grid_F,N_grid);
% % IF = grid_SSM_interp1_nc(ICTT,Grid_F,N_grid);
% % 
% % RTF = @(xa) convert_fit_to_value_x(RF,xa,R0);
% % WTF = @(xa) convert_fit_to_value_x(WF,xa,W0);
% % 
% % ICTTF = @(xa) convert_fit_to_value_x_na(IF,xa,zeros(n,1));
% % 
% % Load_Check = cell2mat(ICTT);
% % [dXdF, Xtrunc,F_trunc] = ftd(Load_Check, Grid_F);
% % Boad_Check = cell(1,max(size(F_trunc)));
% % for ind = 1:max(size(F_trunc))
% % Boad_Check{1,ind} = dXdF(:,ind);
% % end    
% % DIF = grid_SSM_interp1_nc(Boad_Check,F_trunc,max(size(F_trunc)));
% % DICTTF = @(xa) convert_fit_to_value_x_na(DIF,xa,zeros(n,1));
% % % % 
% % Omega = 0.01;
% % q0 = 0.6*exp(1i*0.5);
% % q0 = [q0;conj(q0)];
% % target = @(alpha) {Force_Lorenz(alpha)};
% % z0 = reduced_to_full_traj_adiabatic_sgrid(0,q0,WTF,Omega,target);
% % tf = 60/Omega;
% % nsteps = 2000;
% % outdof = [29];
% % 
% % tic
% % traj3 = transient_traj_on_auto_ssm_adiabatic_sgrid(DS, resonant_modes, WTF, RTF, Omega,target,tf, nsteps, outdof, [], q0);
% % toc
% % % % % 
% 
% % % % % % % 
% % tic
% % ctspan = linspace(0,tf,2001);
% % IC = (z0 +[ICTTF(target(0));DICTTF(target(0))*Omega*DForce_Lorenz(0)]);
% % [tSP,SP_Traj] = ode45(Gh,ctspan,IC);
% % toc
% % 
% % % % 
% % % % tic
% % % % ctspan = linspace(0,tf,2001);
% % % % IC = (z0 +[ICTTF(target(0));DICTTF(target(0))*Omega*DForce_Lorenz(0)])+1;
% % % % [tSP_Off,SP_Traj_Off] = ode45(Gh,ctspan,IC);
% % % % toc
% % % % 
% % 
% %  shift_traj = shift_to_target(DICTTF,ICTTF,traj3,Omega,target,DForce_Lorenz);
% % % % % 
% % traj_actual = traj3.phy;
% % for ind = 1:max(size(traj3.time))
% %    traj_actual(ind,:) =[ICTTF(target(Omega*traj3.time(ind)));DICTTF(target(Omega*traj3.time(ind)))*Omega*DForce_Lorenz(Omega*traj3.time(ind))].';
% % end
% % % % 
% % figure 
% % plot(tSP,SP_Traj(:,1),'-','LineWidth',3,'color','black');
% % hold on
% % plot(traj3.time,shift_traj(:,1),'--','LineWidth',3,'color','red');
% % hold on 
% % plot(traj3.time,traj_actual(:,1),'o','color','blue');
% % xlabel('$t \,[$s$]$','Interpreter','latex');
% % ylabel('$q_1 \,[$m$]$','Interpreter','latex');
% 
% 
% % % % 
% % % 
% % % figure 
% % % plot(linspace(0,60/(0.1),3001),xi0,'-*','LineWidth',3,'color','black');
% % % hold on
% % % plot(linspace(0,60/(0.05),3001),xi1,'-o','LineWidth',3,'color','red');
% % % hold on 
% % % plot(linspace(0,60/(0.01),3001),xi2,'-square','color','blue');
% % % hold on
% % % plot(linspace(0,60/(0.001),3001),xi3,'-+','color','green');
% % % 
% % % xlabel('$t \,[$s$]$','Interpreter','latex');
% % % ylabel('$|\xi|$','Interpreter','latex');
% % % 
% % % axis([0 600 0 0.9])
% % % 
% % % 
% % % figure 
% % % plot(linspace(0,60,3001),xi0,'-*','LineWidth',3,'color','black');
% % % hold on
% % % plot(linspace(0,60,3001),xi1,'-o','LineWidth',3,'color','red');
% % % hold on 
% % % plot(linspace(0,60,3001),xi2,'-square','color','blue');
% % % hold on
% % % plot(linspace(0,60,3001),xi3,'-+','color','green');
% % % 
% % % xlabel('$\alpha$','Interpreter','latex');
% % % ylabel('$|\xi|$','Interpreter','latex');
% % % 
% % % axis([0 60 0 0.9])
% 
% % % % 
% % % % tic
% % % % tspan = linspace(0,10,1000);
% % % % IC = (z0*0 +[ICTTF(target(0));0*DICTTF(target(0))*Omega*DForce_Lorenz(0)]);
% % % % [tBeam,Beam_SSM_Traj] = ode15s(Gh,tspan,IC);
% % % % toc
% % % % 
% % % % tic
% % % % tspan = linspace(0,tf,nsteps);
% % % % Off_epsilon = 0.02;
% % % % IC = IC + Off_epsilon;
% % % % [Off_tBeam,Off_Beam_SSM_Traj] = ode45(Gh,tspan,IC);
% % % % toc
% % % % 
% 
% 
% 
% 
% 
% 
% 
