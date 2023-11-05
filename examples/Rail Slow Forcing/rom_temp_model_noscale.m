function dy = rom_temp_model_noscale(t,zk,l1,gamma,epsilon,xi_sol1)

z = zk(1,1)+1i*zk(2,1);

% alpha = 2;
% beta = 0.1;
% a = 1;

% T = pi/3* sin(Omega*t);
% ls = -1/2*(a*beta+1i*sqrt(a)*sqrt(4-a*beta^2));
% lf = -1/2*((a+2*b)*beta+1i*sqrt(a+2*b)*sqrt(4-(a+2*b)*beta^2));

% es = 1/(sqrt(2*(1+a)*(1+sin(alpha*T))))*[(1+sin(alpha*T))/(sqrt(a))*conj(ls),cos(alpha*T)/sqrt(a)*conj(ls),sqrt(a)*(1+sin(alpha*T)),sqrt(a)*cos(alpha*T)];
% ef = 1/(sqrt(2*(1+a+2*b)*(1-sin(alpha*T))))*[(sin(alpha*T)-1)/(sqrt(a+2*b))*conj(lf),cos(alpha*T)/(sqrt(a+2*b))*conj(lf),sqrt(a+2*b)*(sin(alpha*T)-1),sqrt(a+2*b)*cos(alpha*T)];
% V = [es.',conj(es.'),ef.',conj(ef.')];
% VI = inv(V);
% K1 = VI(1,:)*[0;0;0;-1i/2];
% phaseK = K1/(1i*abs(K1));
% V = phaseK*V;
% VI = inv(V);
% K1 = VI(1,:)*[0;0;0;-1i/2];

% normalise V and VI such that K is purely imaginary
% lf1 = lf;
% ls1 = ls;
% 
% % FPT = VI*[0;0;fp(t);fp(t)];
% 
% g =-gamma;
% gamma_tilde = 0.4;
% g_t = -gamma_tilde;
% 
% 
% ap = V(1,1);
% bp = V(1,3);
% L1 = VI(1,3);
% L3 = VI(3,3);
% 
% e6 = -(ap^3*g*L3)/(lf1 - 3*ls1);
% e7 = -(conj(ap)^3*g*L3)/(lf1 - 3*conj(ls1)); 
% e8 = (3* ap^2 *g*L3*conj(ap))/(-lf1 + 2*ls1 + conj(ls1)); 
% e9 = (3*ap*g*L3*conj(ap)^2)/(-lf1 + ls1 + 2*conj(ls1));

% e15 = (- (3 *ap^6 *g^2 *L1* L3)/(lf1 - 3 *ls1) + (3 *ap^5 *bp *g^2 *L3^2)/(lf1 - 3 *ls1)  + (3 *ap^5 *g *L3 *conj(bp)* g* conj(L3))/(-3 *ls1 + conj(lf1)) + (3 *ap^5 *g^2* L3* conj(ap)* conj(L1))/(-lf1 + 2*ls1 + conj(ls1)) )/(lf1  - 5*ls1); 
% 
% e16  = ((3*ap^3*g^2*L3^2*conj(ap)^2*conj(bp))/(lf1 - 3*ls1)  - (3*g^2*L3*conj(ap)^6*conj(L1))/(lf1 - 3*conj(ls1)) + 3*g^2*L3*conj(ap)^5*((bp* L3)/(lf1 - 3*conj(ls1)) + (ap* L1)/(-lf1 + ls1 + 2*conj(ls1))))/(lf1 - 5*conj(ls1));
% 
% e17 = (- (9*ap^4*g*L3*conj(bp)*conj(ap*g*L3))/(2*ls1 - conj(lf1) + conj(ls1)) + 3*ap^4*g^2*L3*conj(ap)^2* conj(L1)*(3/(-lf1 + 2*ls1 + conj(ls1)) + 2/(-lf1 + ls1 + 2*conj(ls1))) + 3*ap^4*g*L3*conj(ap)*((2*conj(bp)*conj(g*L3))/(-3*ls1 + conj(lf1)) + (g*(-5*ap*L1*lf1 + 5*bp*L3*lf1 + 12*ap*L1*ls1 - 13*bp*L3*ls1 + (3*ap*L1 - 2*bp*L3)*conj(ls1)))/((lf1 - 3*ls1)*(lf1 - 2*ls1 - conj(ls1)))))/(lf1 - 4*ls1 -conj(ls1));
% 
% 
% 
% e18 = ( (6*ap^4*g^2*L3^2*conj(ap)*conj(bp))/(lf1 - 3*ls1) + 3*ap*g^2*L3*conj(ap)^5*conj(L1)*(-(3/(lf1 - 3*conj(ls1))) + 2/(-lf1 + ls1 + 2*conj(ls1))) + 3*ap*g*L3*conj(ap)^4*(-((3*conj(bp)*conj(g*L3))/(ls1 - conj(lf1) + 2*conj(ls1))) + g *((2* bp* L3)/(lf1 - 3*conj(ls1)) + (3*bp*L3)/(lf1 - ls1 - 2*conj(ls1)) + (2*ap*L1)/(-lf1 + 2*ls1 + conj(ls1)) + (3*ap*L1)/(-lf1 + ls1 + 2*conj(ls1)))))/(lf1 - ls1 - 4*conj(ls1));
% 
% 
% e19 =(-((18*ap^3*g*L3*conj(ap)*conj(bp)*conj(ap*g*L3))/(2*ls1 - conj(lf1) + conj(ls1)))  + 3*ap^3*g^2*L3*conj(ap)^3*conj(L1)*(3/(-lf1 + 2*ls1 + conj(ls1)) + 6/(-lf1 + ls1 + 2*conj(ls1)) + 1/(-lf1 + 3*conj(ls1))) + 3*ap^3*g*L3*conj(ap)^2*(conj(bp)*conj(g*L3)*(1/(-3*ls1 + conj(lf1)) - 3/(ls1 - conj(lf1) + 2*conj(ls1))) + g*(bp*L3*(1/(lf1 - 3*ls1) + 3/(lf1 - ls1 - 2*conj(ls1)) + 6/(lf1 - 2*ls1 - conj(ls1))) + ap*L1* (-(3/(lf1 - 3* ls1)) + 6/(-lf1 + 2 *ls1 + conj(ls1)) + 1/(-lf1 + ls1 + 2*conj(ls1))))))/(lf1 - 3*ls1 + 2*conj(ls1));
% 
% 
% e20 = ((3*ap^5*g^2*L3^2*conj(bp))/(lf1 - 3*ls1) - (9*ap^2*g*L3*conj(ap)^2* conj(bp) *conj(ap*g*L3))/(2*ls1 - conj(lf1) + conj(ls1)) + 3*ap^2 *g^2* L3*conj(ap)^4*conj(L1)*(-(3/(lf1 - 3*conj(ls1))) +  1/(-lf1 + 2*ls1 + conj(ls1)) + 6/(-lf1 + ls1 + 2*conj(ls1))) + 3*ap^2* g *L3* conj(ap)^3*(bp*g*L3*(1/(lf1 - 3 *conj(ls1)) + 6/(lf1 - ls1 - 2*conj(ls1)) + 3/( lf1 - 2 *ls1 - conj(ls1))) - (6 *conj(bp)*conj(g*L3))/( ls1 - conj(lf1) + 2*conj(ls1)) + ap*g* L1* (1/(-lf1 + 3* ls1) + 6/(-lf1 + 2* ls1 + conj(ls1)) + 3/(-lf1 + ls1 + 2 *conj(ls1)))))/(lf1 -2*ls1 - 3*conj(ls1));

% p = epsilon^(3/M)*(e6*z^3+e7*conj(z)^3+e8*z^2*conj(z)+e9*conj(z)^2*z);
% 
% dz = ls*z+epsilon^(-1/M)*L1*g*(epsilon^(3/M)*(ap*z+conj(ap*z))^3+epsilon^(2/M)*3*(ap*z+conj(ap*z))^2*(bp*p+conj(bp*p)))+ epsilon^(1-1/M)*FPT(1);
% dz = ls*z+L1*g*((ap*z+conj(ap*z))^3+3*(ap*z+conj(ap*z))^2*(bp*p+conj(bp*p)))+epk*VI(1,4)*sin(Opema*t);

% p = (e6*z^3+e7*conj(z)^3+e8*z^2*conj(z)+e9*conj(z)^2*z);
% p1 = h021(t)*conj(z)^2 +h201(t)*(z)^2 + h111(t)*z*conj(z); 
% p2 = h012(t)*conj(z) + h102(t)*z;
% sol1 = XI_hyperbolic_solution{1,1}(t);
% sol3 = XI_hyperbolic_solution{1,3}(t);
% sol11 = sol1(1);
% sol13 = sol3(1);
% OP1 = epsilon*(bp*p1+conj(bp*p1));
% dz = ls*z+L1*g*((ap*z+conj(ap*z))^3+3*(ap*z+conj(ap*z))^2*(bp*p+conj(bp*p)))+ L1*g*3*epsilon*sol11*(ap*z+conj(ap*z))^2+L1*g*3*epsilon^3*sol13*(ap*z+conj(ap*z))^2+L1*g*3*epsilon*sol11*2*(ap*z+conj(ap*z))*OP1+L1*g*3*epsilon^2*sol11^2*(ap*z+conj(ap*z)+bp*p+conj(bp*p)); % + L1*g_t*((ap*z+conj(ap*z))^2 + 2*(ap*z+conj(ap*z))*(bp*p+conj(bp*p)));
% dz = dz + 3*L1*g*(ap*z+conj(ap*z))^2*epsilon^2*(bp*p2+conj(bp*p2));

% O1 = bp*p+conj(bp*p) + epsilon*(bp*p1+conj(bp*p1)) + epsilon^2*(bp*p2+conj(bp*p2));
% dz = ls*z + L1*g*((ap*z+conj(ap*z))^3+3*(ap*z+conj(ap*z))^2*O1 + 3*epsilon*sol11*(ap*z+conj(ap*z))*2*O1 + 3*epsilon*sol11*(ap*z+conj(ap*z))^2 + 3*epsilon^3*sol13*(ap*z+conj(ap*z))^2 + 3*epsilon^2*sol11^2*(ap*z+conj(ap*z)+O1) );

Np = 1;

u1 = z;
u2 = conj(z);
xi1 = xi_sol1(t);
du = (-0.000331578 + 0.00105703*1i)*u1^3*gamma - (0.000500717 - 0.00328552*1i)*u1^2*u2*gamma + (4.95865*10^-6 + 0.00332345*1i)*u1*u2^2*gamma  + (0.000170173 + 0.00109467*1i)*u2^3*gamma + u1*l1 - (0.00258524 - 0.0111618*1i)* u1^2*gamma*epsilon*xi1 - (0.00171407 - 0.0228504*1i)* u1* u2*gamma*epsilon*xi1 + (0.000891124 + 0.0114226*1i)* u2^2*gamma*epsilon*xi1 - (0.00595079 - 0.0390469*1i)* u1*gamma*epsilon^2*xi1^2 + (0.0000589313 + 0.0394977*1i)* u2*gamma*epsilon^2*xi1^2;


dy(1:Np,1) = real(du);
dy(Np+1:2*Np,1) = imag(du);


end