clear all
Order = 3;

XI_hyperbolic_solution = cell(1,Order);

c = 0.3;
k = 1;
a = 0.3;
g = 9.8;
m = 1;
M = 4;
cf = 0.3;
scale = 5*a^3;%2*a^3;

K = [-2*a^2*(1 + m/M)*g*4/scale - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];


B = eye(4);
signe = -1;
% quadratic part
F2 = sptensor([4 4 4]);
F3 = sptensor([4 4 4 4]);
F3(3,1,1,1) = -4*g*(1+m/M)/scale;
F3(3,1,3,3) = -16*(1+m/M)*a^4/scale^2;
F2(3,1,1) =signe*-4*g*(1+m/M)*3*a/scale;
F2(3,3,3) = signe*-16*(1+m/M)*a^5/scale^2;

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

%Lorenz force 
%%

epsilon = 0.1;
L = 0;
N_Max = 10000;
tspan = linspace(L,500,N_Max);

loren = @(t,y) lorenz_3D(t,y);
IC = [0;0.3;0.5];
[tk,yy] = ode45(loren,tspan,IC);

sigma = 10;
rho = 28;
beta = 8/3;


F_a_new = yy(:,1)/max(abs(yy(:,1)));

Force_Lorenz = griddedInterpolant(tspan(:),F_a_new(:),'spline');



xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));Force_Lorenz(ta)])));

positive_L = 1;

tcalc = tspan(positive_L:end);
N_Max = max(size(tcalc));
LT = [];
start_time = [1];
xi_sol = cell(1,max(size(start_time)));

tic
sol_1 = [];
parfor i = 1:max(size(tcalc))
Lo = positive_L - 1 + i ;
T = tspan(start_time(1):Lo);
ti = tcalc(i);

INx = xs(ti,T);

Q1x = trapz(T,INx,2);
sol_1  =[sol_1,Q1x];
end    
xi_sol{1,1} = sol_1;
LT = [LT,tspan(start_time(1))];
toc
 start_time=positive_L;
xi_solg=sol_1.'; 
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval = @(xa) xi_full(xa);
save('Solution_1.mat','xi_full_eval')
XI_hyperbolic_solution{1,1} = xi_full_eval;
xi_sol0=xi_sol{1,1};
xi_sol11 = griddedInterpolant(tcalc(:),(xi_sol0(1,:)).','spline');
% xi_sol2 = @(xq) interp1(tcalc(:),(xi_sol0(2,:)).',xq,'spline');
xi_sol13 = griddedInterpolant(tcalc(:),(xi_sol0(3,:)).','spline');
% xi_sol4 = @(xq) interp1(tcalc(:),(xi_sol0(4,:)).',xq,'spline');
% xi_net_sol = @(t) [xi_sol1(t);xi_sol2(t);xi_sol3(t);xi_sol4(t)];
save('Xi_Sol13_stable.mat','xi_sol11','xi_sol13')


xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));signe*(-4*g*(1+m/M)*3*a*xi_sol11(ta).^2/scale -16*(1+m/M)*a^5*xi_sol13(ta).^2/scale^2);0*ones(1,max(size(ta)))])));

% positive_Lvec = find(tspan<-2999 & tspan>-3000);
% positive_L = positive_Lvec(1)-1;
positive_L = 1;

tcalc = tspan(positive_L:end);
N_Max = max(size(tcalc));
LT = [];

xi_sol = cell(1,max(size(start_time)));

tic
sol_1 = [];
parfor i = 1:max(size(tcalc))
Lo = positive_L - 1 + i ;
T = tspan(start_time:Lo);
ti = tcalc(i);

INx = xs(ti,T);

Q1x = trapz(T,INx,2);
sol_1  =[sol_1,Q1x];
end    
xi_sol{1,1} = sol_1;
LT = [LT,tspan(start_time(1))];
toc
xi_sol{1,1} = sol_1;

xi_solg=sol_1.'; 
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval2 = @(xa) xi_full(xa);
save('Solution_2Stable.mat','xi_full_eval2')
XI_hyperbolic_solution{1,2} = xi_full_eval2;

xi_sol0=xi_sol{1,1};
xi_sol21 = griddedInterpolant(tcalc(:),(xi_sol0(1,:)).','spline');
xi_sol23 = griddedInterpolant(tcalc(:),(xi_sol0(3,:)).','spline');





xs=@(t,ta) real(V*([exp(D(1)*(t-ta));conj(exp(D(1)*(t-ta)));exp(D(3)*(t-ta));conj(exp(D(3)*(t-ta)))].*(VI*[0*ones(1,max(size(ta)));0*ones(1,max(size(ta)));-4*g*(1+m/M)*xi_sol11(ta).^3/scale-16*(1+m/M)*a^4*xi_sol13(ta).^2.*xi_sol11(ta)/scale^2+signe*(-4*g*(1+m/M)*3*a*2*xi_sol11(ta).*xi_sol21(ta)/scale -16*(1+m/M)*a^5*2*xi_sol13(ta).*xi_sol23(ta)/scale^2);0*ones(1,max(size(ta)))])));

positive_L = 1;

tcalc = tspan(positive_L:end);
N_Max = max(size(tcalc));
LT = [];

xi_sol = cell(1,max(size(start_time)));

tic
sol_1 = [];
parfor i = 1:max(size(tcalc))
Lo = positive_L - 1 + i ;
T = tspan(start_time:Lo);
ti = tcalc(i);

INx = xs(ti,T);

Q1x = trapz(T,INx,2);
sol_1  =[sol_1,Q1x];
end    
xi_sol{1,1} = sol_1;
LT = [LT,tspan(start_time(1))];
toc

xi_solg=sol_1.'; 
xi_full = griddedInterpolant(tcalc(:),xi_solg,'spline');
xi_full_eval3 = @(xa) xi_full(xa);
save('Solution_3Stable.mat','xi_full_eval3')
XI_hyperbolic_solution{1,3} = xi_full_eval3;
save('Mixed_Mode_Stable_minus_a.mat','XI_hyperbolic_solution')


% 