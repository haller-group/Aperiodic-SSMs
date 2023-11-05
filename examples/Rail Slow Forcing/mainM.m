
% [M,C,K,fnl,f_0,outdof,PlotFieldonDefMesh] = build_model(nElements);
clear all
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
% F3(3,1,3,3) = -16*(1+m/M)*a^4;
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
%Lorenz force 
%%
% epsilon = 0.1;
L = 0;
N_Max = 10000;
tspan = linspace(L,500,N_Max);

% probing_window = find(tspan<=200 & tspan>199);
% T_pos = probing_window(end);
% Time_pos = tspan(T_pos);

loren = @(t,y) lorenz_3D(t,y);
IC = [0;0.3;0.5];
[tk,yy] = ode45(loren,tspan,IC);

sigma = 10;
rho = 28;
beta = 8/3;


F_a_new = yy(:,1)/max(abs(yy(:,1)));

% Force_Lorenz = @(xq) interp1(tspan(:),F_a_new(:),xq,'spline');

Force_Lorenz = griddedInterpolant(tspan(:),F_a_new(:),'spline');

% ty = linspace(0,50,100);
% F = [1*ones(1,50),2*ones(1,50)];
% figure 
% plot(ty(1:50),F(1:50),'-','LineWidth',3,'color','blue')
% hold on 
% plot(ty(51:end),F(51:end),'-','LineWidth',3,'color','blue')
% axis([0 50 0 3])
% xlabel('$t \,[$s$]$','Interpreter','latex');
% ylabel('$F(t) \,[$N$]$','Interpreter','latex');





%%
% order epsilon solution -unstable


% order epsilon solution -unstable
% AP = [0,0,1,0;0,0,0,1;-2*k,k,-c,c;k,-2*k,c,-2*c];

xs_neg=@(t,ta) real(V*([zeros(1,max(size(ta)));(exp(D(2)*(t-ta)));exp(D(3)*(t-ta));(exp(D(4)*(t-ta)))].*(inv(V)*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta)])));
xs_pos=@(t,ta) real(V*([exp(D(1)*(t-ta));zeros(1,max(size(ta)));zeros(1,max(size(ta)));zeros(1,max(size(ta)))].*(inv(V)*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta)])));
% 
neg_L =1;


% 
% 
tcalc = tspan;
N_Max = max(size(tcalc));
LT = [];
start_time = 1;
xi_sol = cell(1,max(size(start_time)));

tic
sol_1 = [];
parfor i = 1:max(size(tcalc))
Lo_neg = neg_L - 1 + i ;
T_neg = tspan(start_time(1):Lo_neg);
T_pos = tspan(Lo_neg:end);

ti = tcalc(i);

INx_neg = xs_neg(ti,T_neg);
INx_pos = xs_pos(ti,T_pos);

Q1x_neg = trapz(T_neg,INx_neg,2);
Q1x_pos = trapz(T_pos,INx_pos,2);

sol_1  =[sol_1,Q1x_neg-Q1x_pos];
end    
xi_sol{1,1} = sol_1;
LT = [LT,tspan(start_time(1))];
toc

xi_solg=sol_1.'; 
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval = @(xa) xi_full(xa);
save('Order_Ep_1_mixed_unstable.mat','xi_full_eval')

xi_sol0=xi_sol{1,1};
xi_sol1 = griddedInterpolant(tcalc(:),(xi_sol0(1,:)).','spline');
xi_sol2 = @(xq) interp1(tcalc(:),(xi_sol0(2,:)).',xq,'spline');
xi_sol3 = griddedInterpolant(tcalc(:),(xi_sol0(3,:)).','spline');
xi_sol4 = @(xq) interp1(tcalc(:),(xi_sol0(4,:)).',xq,'spline');
xi_net_sol = @(t) [xi_sol1(t);xi_sol2(t);xi_sol3(t);xi_sol4(t)];
save('Order_Ep_1_mixed_unstable13.mat','xi_sol1','xi_sol3')

load('Order_Ep_1_mixed_unstable.mat','xi_full_eval')
load('Order_Ep_1_mixed_unstable13.mat','xi_sol1','xi_sol3')

start_time = neg_L;


xs_neg=@(t,ta) real(V*([zeros(1,max(size(ta)));(exp(D(2)*(t-ta)));exp(D(3)*(t-ta));(exp(D(4)*(t-ta)))].*(inv(V)*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));(-4*g*(1+m/M)*xi_sol1(ta).^3-16*(1+m/M)*a^4.*xi_sol1(ta).*xi_sol3(ta).^2);0*ones(1,max(size(ta)))])));
xs_pos=@(t,ta) real(V*([exp(D(1)*(t-ta));zeros(1,max(size(ta)));zeros(1,max(size(ta)));zeros(1,max(size(ta)))].*(inv(V)*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));(-4*g*(1+m/M)*xi_sol1(ta).^3-16*(1+m/M)*a^4.*xi_sol1(ta).*xi_sol3(ta).^2);0*ones(1,max(size(ta)))])));




N_Max = max(size(tcalc));
LT = [];

xi3 = cell(1,max(size(start_time)));

tic
sol_1 = [];
parfor i = 1:max(size(tcalc))
Lo_neg = neg_L - 1 + i ;
T_neg = tspan(start_time:Lo_neg);
T_pos = tspan(Lo_neg:end);

ti = tcalc(i);

INx_neg = xs_neg(ti,T_neg);
INx_pos = xs_pos(ti,T_pos);

Q1x_neg = trapz(T_neg,INx_neg,2);
Q1x_pos = trapz(T_pos,INx_pos,2);

sol_1  =[sol_1,Q1x_neg-Q1x_pos];
end    
xi3{1,1} = sol_1;
LT = [LT,tspan(start_time(1))];
toc

xi_solg=sol_1.'; 
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval3 = @(xa) xi_full(xa);
save('Order_Ep_3_mixed_unstable.mat','xi_full_eval3')

%%
%Compute coefficients of SSM and reduced order models

% l1 = -3.6932900912556486;
% l2 = 3.3037897062314343;
% l3 = -0.03024980748789345 + 0.44707282565262196*1i;
% % a = 0.5
% h300 = ((-0.003309525166448172 + 0.02244712920639276*1i))/(l3-3*l1);
% h030 = ((0.0037438910548629844 - 0.025393251906695717*1i))/(l3-3*l2);
% h120 = ((0.0015513520137821193 - 0.010522173830555212*1i))/(l3-l1-2*l2);
% h210 = ((-0.000778204351452102 + 0.005278235622172289*1i))/(l3-2*l1-l2);

l3 = -0.030061132181587 + 0.4464257565641846*1i;
l1 = -7.194163190718024;
l2 = 6.804285455081199;
% a = 1;
h300 = ((-0.00038667273282006853 + 0.0027979839293201573*1i))/(l3-3*l1);
h030 = ((0.00040962950175648337 - 0.0029641002988007393*1i))/(l3-3*l2);
h120 = ((0.00035210505940452244 - 0.002547850453433994*1i))/(l3-l1-2*l2);
h210 = ((-0.000292279522327395 + 0.002114949767409561*1i))/(l3-2*l1-l2);

hc300 = conj(h300);
hc030 = conj(h030);
hc120 = conj(h120);
hc210 = conj(h210);
h500 =(1.389965047288079*h210 + (4.4040673648363295 - 0.00010373937608599797*1i)*h300 + (0.000011284585006269063 - 0.00010459931933348712*1i)*hc300)/(l3-5*l1);

h050 = ((-4.4174221584054845 - 0.0000971855068212439*1i)*h030 - 1.5551727992695028*h120 - (0.000011904761269644814 + 0.00010432442367424514*1i)*hc030 )/(l3-5*l2);

h320 = (4.169895141864237*h030 + (3.569337298774463 - 0.00010373937608599797*1i)*h120 + (0.9536361398675327 - 0.0001483720984409183*1i)*h210 - ...
  (4.010297943983738 + 0.00009718550682124392*1i)*h300 + (0.000011284585006269063 - 0.00010459931933348715*1i)*hc120 - (6.460597498676684*10^-7 + 0.00015433182163900706*1i)*hc210 - ...
  (0.000011904761269644819 + 0.00010432442367424517*1i)*hc300)/(l3-3*l1 -2*l2) ;

h230 = ((3.1519722657435296 - 0.00010373937608599797*1i)*h030 - (1.4217186338734653 + 0.00014837209844091828*1i)*h120 - (4.146006015457653 + 0.00009718550682124387*1i)*h210 - ...
  4.665518397808508*h300 + (0.000011284585006269056 - 0.00010459931933348715*1i)*hc030 - (6.460597498676684*10^-7 + 0.00015433182163900698*1i)*hc120 - ...
  (0.000011904761269644819 + 0.00010432442367424518*1i)*hc210)/(l3-3*l2-2*l1);

h410 = (2.779930094576158*h120 + (3.986702331805396 - 0.00010373937608599797*1i)*h210 + (3.3289909136085303 - 0.0001483720984409183*1i)*h300 + ...
  (0.000011284585006269063 - 0.00010459931933348712*1i)*hc210 - (6.460597498676684*10^-7 + 0.000154331821639007*1i)*hc300)/(l3-4*l1-l2);

h140 = ((-3.7970734076144628 - 0.0001483720984409183*1i)*h030 - (4.281714086931569 + 0.00009718550682124391*1i)*h120 - 3.110345598539005*h210 - ...
  (6.460597498676683*10^-7 + 0.000154331821639007*1i)*hc030 - (0.000011904761269644817 + 0.00010432442367424514*1i)*hc120)/(l3-4*l2-l1);


h300 = ((-0.000017476819444004 + 0.00012646317102144345*1i))/(l3-3*l1);
h030 = ((0.00002058722597287155 - 0.00014897023382349944*1i))/(l3-3*l2);
h120 = ((-0.00005847998809831887 + 0.0004231642238969843*1i))/(l3-l1-2*l2);
h210 = ((0.00005537266996673961 - 0.00040067950889772515*1i))/(l3-2*l1-l2);

% 
% 
h_coeff_pos=@(t,ta) [exp((l3-l2-l1)*(t-ta)).*-1.*((0.0045528085943358785 - 0.03294435888281321*1i).* xi_sol1(ta) + (0.00004266803042538553 - 0.00030874807891231816*1i).* ...
xi_sol3(ta)) ;
exp((l3-2*l1)*(t-ta)).*-1.*((0.003062413331418653 - 0.022159781538649215*1i).*xi_sol1(ta) - ...
(0.0007454913958689365 - 0.0053944143665758835*1i).*xi_sol3(ta)) ;
 exp((l3-l1)*(t-ta)) ...
.*-1.*(- (0.0027660226580715234 - 0.020015083269450215*1i).* xi_sol1(ta).^2 + (0.005414753303779644 - 0.0391814354602*1i).* xi_sol1(ta).* xi_sol3(ta) - (0.0003763296133430644 - 0.00272314058087758*1i).*xi_sol3(ta).^2) ;
];

% h111 exp((l3-l2-l1)*(t-ta)).*(0.0008043804376629333 - 0.00582053852420905*1i).*xi_sol1(ta)
% h201 -1.*exp((l3-2*l1)*(t-ta)).*(0.0003808199518400002 - 0.00275563290289885*1i).*xi_sol1(ta)
% h102 exp((l3-l1)*(t-ta)).*(0.002766022658071523 - 0.02001508326945021*1i).*xi_sol1(ta).^2

% h021 -1.*exp((l3-2*l2)*(t-ta)).*(0.00042475970952190155 - 0.0030735832661311887*1i).*xi_sol1(ta)
% h012 -1.*exp((l3-l2)*(t-ta)).*(0.0029212420535412995 - 0.02113825885743593*1i).*xi_sol1(ta).^2

h_coeff_neg=@(t,ta) [exp((l3-2*l2)*(t-ta)).*-1.*((0.0031003587156133655 - 0.022434356304742407*1i).*xi_sol1(ta) + (0.0007864452553481337 - 0.005690758615705421*1i).* ...
xi_sol3(ta)) ;
exp((l3-l2)*(t-ta)) ...
.*-1.*(+(0.0029212420535412995- 0.02113825885743593*1i).*...
xi_sol1(ta).^2 + (0.005408697909029281- 0.03913761831003883*1i).* xi_sol1(ta) .*xi_sol3(ta) + (0.0003974478984409931 - 0.0028759535860457047*1i ).* xi_sol3(ta).^2)
];

N_Max = max(size(tcalc));
sol_1r_neg = [];
sol_1i_neg = [];

sol_1r_pos = [];
sol_1i_pos = [];

tic
parfor i = 1:max(size(tcalc))
Lo = neg_L - 1 + i ;
T_neg = tspan(1:Lo);
T_pos = tspan(Lo:end);
ti = tcalc(i);

INx_neg = h_coeff_neg(ti,T_neg);
INxr_neg = real(INx_neg);
INxi_neg = imag(INx_neg);
Q1xr_neg = trapz(T_neg,INxr_neg,2);
Q1xi_neg = trapz(T_neg,INxi_neg,2);
sol_1r_neg  =[sol_1r_neg,Q1xr_neg];
sol_1i_neg  =[sol_1i_neg,Q1xi_neg];

INx_pos = h_coeff_pos(ti,T_pos);
INxr_pos = real(INx_pos);
INxi_pos = imag(INx_pos);
Q1xr_pos = trapz(T_pos,INxr_pos,2);
Q1xi_pos = trapz(T_pos,INxi_pos,2);
sol_1r_pos  =[sol_1r_pos,-1*Q1xr_pos];
sol_1i_pos  =[sol_1i_pos,-1*Q1xi_pos];

end    
toc

h111r = griddedInterpolant(tcalc(:),(sol_1r_pos(1,:)).','spline');
h111i = griddedInterpolant(tcalc(:),(sol_1i_pos(1,:)).','spline');
h111 = @(xa) h111r(xa) + 1i*h111i(xa);

h201r = griddedInterpolant(tcalc(:),(sol_1r_pos(2,:)).','spline');
h201i = griddedInterpolant(tcalc(:),(sol_1i_pos(2,:)).','spline');
h201 = @(xa) h201r(xa) + 1i*h201i(xa);

h102r = griddedInterpolant(tcalc(:),(sol_1r_pos(3,:)).','spline');
h102i = griddedInterpolant(tcalc(:),(sol_1i_pos(3,:)).','spline');
h102 = @(xa) h102r(xa) + 1i*h102i(xa);

h021r = griddedInterpolant(tcalc(:),(sol_1r_neg(1,:)).','spline');
h021i = griddedInterpolant(tcalc(:),(sol_1i_neg(1,:)).','spline');
h021 = @(xa) h021r(xa) + 1i*h021i(xa);

h012r = griddedInterpolant(tcalc(:),(sol_1r_neg(2,:)).','spline');
h012i = griddedInterpolant(tcalc(:),(sol_1i_neg(2,:)).','spline');
h012 = @(xa) h012r(xa) + 1i*h012i(xa);

%%

save('Mixed_Mode_Coeff_SSM.mat','h030','h300','h120','h210','h111','h201','h021','h012','h102')
load('Mixed_Mode_Coeff_SSM.mat','h030','h300','h120','h210','h111','h201','h021','h012','h102')

%%
% 
load('Order_Ep_1_mixed_unstable.mat','xi_full_eval')
load('Order_Ep_1_mixed_unstable13.mat','xi_sol1')
% 
load('Order_Ep_3_mixed_unstable.mat','xi_full_eval3')
% 
Ma1 = load('Mixed_Mode_Stable_plus_a.mat','XI_hyperbolic_solution');
XI_hyperbolic_solution_a1 = Ma1.('XI_hyperbolic_solution');

Ma1 = load('Mixed_Mode_Stable_minus_a.mat','XI_hyperbolic_solution');
XI_hyperbolic_solution_a2 = Ma1.('XI_hyperbolic_solution');
% 
XI_hyperbolic_solution = cell(1,3);
XI_hyperbolic_solution{1,1} = @(xa) xi_full_eval(xa);
XI_hyperbolic_solution{1,2} = {};
XI_hyperbolic_solution{1,3} = @(xa) xi_full_eval3(xa);
% 
% 
% 
epsilon = 0;
ctspan = linspace(0,500,8000);
Unstable_Solution = epsilon*xi_full_eval(ctspan) + epsilon^3*xi_full_eval3(ctspan);

Stable_solution1 = epsilon*XI_hyperbolic_solution_a1{1,1}(ctspan) + epsilon^2*XI_hyperbolic_solution_a1{1,2}(ctspan) + epsilon^3*XI_hyperbolic_solution_a1{1,3}(ctspan);
Stable_solution2 = epsilon*XI_hyperbolic_solution_a2{1,1}(ctspan) + epsilon^2*XI_hyperbolic_solution_a2{1,2}(ctspan) + epsilon^3*XI_hyperbolic_solution_a2{1,3}(ctspan);
Stable_solution1(:,1) = Stable_solution1(:,1)+a;
Stable_solution1(:,2) = Stable_solution1(:,2)+a*m/(M+m);

Stable_solution2(:,1) = Stable_solution2(:,1)-a;
Stable_solution2(:,2) = Stable_solution2(:,2)-a*m/(M+m);

% 
proj = VI_new*(Unstable_Solution(1,:).');
u0 = [proj(1);proj(2)+0.1];
[y0,modal0] = compute_SSM_phy(XI_hyperbolic_solution,0,V_new,VI_new,u0,h030,h300,h120,h210,h111,h201,h021,h012,h102,h500,h050,h410,h140,h230,h320,epsilon);

RomH = @(t,uk) rom_temp_model_noscale(t,uk,l1,l2,V_new,VI_new,epsilon,xi_sol1,xi_sol3,h300,h030,h120,h210);
tic
IC = u0;
[tred,ROM_SOL] = ode45(RomH,ctspan,IC);
toc

[y_sol,modal_sol] = compute_SSM_phy(XI_hyperbolic_solution,tred.',V_new,VI_new,ROM_SOL.',h030,h300,h120,h210,h111,h201,h021,h012,h102,h500,h050,h410,h140,h230,h320,epsilon);

GH = @(t,y) Full_Order_Model_PM(t,y,Force_Lorenz,epsilon);
tic
IC = y0;
[tSP,Full_Traj] = ode45(GH,ctspan,IC);
toc
NMTE = sum(sqrt(sum((y_sol.' - Full_Traj).^2,2)))/(sqrt(max(sum((Full_Traj).^2,2)))*max(size(ctspan)))


figure
indexR=2;
plot(ctspan,Full_Traj(:,indexR),'-','LineWidth',3,'color',[0,0,0]);
hold on 
plot(ctspan,y_sol(indexR,:),'--','LineWidth',3,'color',[0,1,0]);
hold on
plot(ctspan,Unstable_Solution(:,indexR),'-','LineWidth',3,'color',[1,0,0,0.3]);
hold on 
plot(ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.6]);
hold on 
plot(ctspan,Stable_solution2(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.6]);
hold on
legend('Full Order Model Response','Reduced Order Model','O(\epsilon^3) Unstable Solution','O(\epsilon^3) Stable 1','O(\epsilon^3) Stable 2')
title('\epsilon = 5')
ylabel('$X_{R}$ (m)','interpreter' , 'latex')
xlabel('$t$ (s)','interpreter' , 'latex')


% axes( 'Position' ,[.3 .6 .6 .1])
% box on
% plot(ctspan,y_sol(indexR,:),'-','LineWidth',3,'color',[0,1,0,1]);
% hold on 
% plot(ctspan,Full_Traj(:,indexR),'--','LineWidth',3,'color',[0,0,0,1]);
% ylabel('X_R (m)')
% xlabel('t (s)')
% axis([40,80,0.98,1.02])

%%
% figure
% indexR=4;
% plot(ctspan,Unstable_Solution(:,indexR),'-','LineWidth',3,'color',[1,0,0,0.3]);
% hold on 
% plot(ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.3]);
% hold on 
% plot(ctspan,Stable_solution2(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.3]);
% hold on
% plot(ctspan,Full_Traj(:,indexR),'--','LineWidth',3,'color',[0,0,0,1]);
% legend('$O(\epsilon^3)$ Unstable Solution','$O(\epsilon^3)$ Stable 1','$O(\epsilon^3)$ Stable 2','Full Order Model Response')
% title('\epsilon = 0.01 and |\epsilon F|_{max} = 0.32')
% ylabel('$\dot{X}_{CM}$ (m/s)','interpreter' , 'latex')
% xlabel('$t$ (s)','interpreter' , 'latex')

% axes( 'Position' ,[.3 .6 .6 .1])
% box on
% plot(ctspan,Stable_solution1(:,indexR),'-','LineWidth',3,'color',[0,0,1,0.3]);
% hold on 
% plot(ctspan,Full_Traj(:,indexR),'--','LineWidth',3,'color',[0,0,0,1]);
% ylabel('X_R (m)')
% xlabel('t (s)')
% axis([40,80,0.98,1.02])



% xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)))].*(inv(V)*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta);Force_Lorenz(ta)])));
% positive_Lvec = find(tspan<-3999 & tspan>-4000);
% positive_L = positive_Lvec(1)-1;
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
% 
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
