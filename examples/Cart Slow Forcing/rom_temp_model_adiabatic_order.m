function dy = rom_temp_model_adiabatic_order(order,t,zk,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_01,xi_11,xi_21,Valpha,ValphaD,V,A,Force_Lorenz,Dalpha,gamma,epsilon)
u1 = zk(1,1)+1i*zk(2,1);
u2 = zk(1,1)-1i*zk(2,1);
m1 = 1;
m2 = 1;
Mf = 2;
kf = 1;
k = 1;
c = 0.3;
cf = 0.3;
gamma = 0.5;
alpha = epsilon*t;
dia = Dalpha(Force_Lorenz(alpha));
l1 = dia(1);
l1c = conj(l1);
l2 = dia(3);
l3 = dia(5);
xi0 = xi_01(alpha);
xi1 = xi_11(alpha,A);
xi2 = xi_21(alpha,A);
P = Valpha(Force_Lorenz(alpha),V);
PI = inv(P);
PD = ValphaD(alpha,V);
TempCell = num2cell(SSM_Coeff_A_1(alpha));
[h101,h200,h011,h110,h020,f101,f200,f011,f110,f020,h300,h030,h120,h210,f300,f030,f120,f210]=deal(TempCell{:});
TempCell = num2cell(conj([h101,h200,h011,h110,h020,f101,f200,f011,f110,f020,h300,h030,h120,h210,f300,f030,f120,f210]));
[hc011,hc020,hc101,hc110,hc200,fc011,fc020,fc101,fc110,fc200,hc030,hc300,hc210,hc120,fc030,fc300,fc210,fc120] = deal(TempCell{:});
% TempCell = num2cell(SSM_Coeff_A_2(alpha));    
% [h012,h102,h201,h021,h111,f012,f102,f201,f021,f111] = deal(TempCell{:});
if order ==3
du =l1.*u1 + u1.^2.*((-3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 4))./m1 - (3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 4))./Mf + (3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 5))./m1) + ...
 u1.*u2.*((-6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 4))./m1 - (6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 4))./Mf + (6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 5))./m1) + ...
 u2.^2.*((-3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 4))./m1 - (3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 4))./Mf + (3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 5))./m1) + ...
 u1.^3.*(-((gamma.*P(1, 1).^3.*PI(1, 4))./m1) - (gamma.*P(1, 1).^3.*PI(1, 4))./Mf - (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - ...
   (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf + ...
   (gamma.*P(1, 1).^3.*PI(1, 5))./m1 + (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + ...
   (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1) + ...
 u2.^3.*(-((gamma.*P(1, 2).^3.*PI(1, 4))./m1) - (gamma.*P(1, 2).^3.*PI(1, 4))./Mf - (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - ...
   (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (gamma.*P(1, 2).^3.*PI(1, 5))./m1 + (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + ...
   (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1) + ...
 u1.*u2.^2.*((-3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 4))./m1 - (3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 4))./Mf - (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 5))./m1 + (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 + (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1) + ...
 u1.^2.*u2.*((-3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 4))./m1 - (3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 4))./Mf - (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 5))./m1 + (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 + (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1) + ...
 epsilon.*u1.*(-(PD(1, 1).*PI(1, 1)) - PD(2, 1).*PI(1, 2) - PD(3, 1).*PI(1, 3) - (6.*gamma.*xi0.*xi1.*P(1, 1).*PI(1, 4))./m1 - (6.*gamma.*xi0.*xi1.*P(1, 1).*PI(1, 4))./Mf - PD(4, 1).*PI(1, 4) + ...
   (6.*gamma.*xi0.*xi1.*P(1, 1).*PI(1, 5))./m1 - PD(5, 1).*PI(1, 5) - PD(6, 1).*PI(1, 6)) + ...
 epsilon.*u2.*(-(PD(1, 2).*PI(1, 1)) - PD(2, 2).*PI(1, 2) - PD(3, 2).*PI(1, 3) - (6.*gamma.*xi0.*xi1.*P(1, 2).*PI(1, 4))./m1 - (6.*gamma.*xi0.*xi1.*P(1, 2).*PI(1, 4))./Mf - PD(4, 2).*PI(1, 4) + ...
   (6.*gamma.*xi0.*xi1.*P(1, 2).*PI(1, 5))./m1 - PD(5, 2).*PI(1, 5) - PD(6, 2).*PI(1, 6)) + ...
 epsilon.^2.*u2.*(-(h011.*PD(1, 3).*PI(1, 1)) - hc011.*PD(1, 4).*PI(1, 1) - f011.*PD(1, 5).*PI(1, 1) - fc011.*PD(1, 6).*PI(1, 1) - h011.*PD(2, 3).*PI(1, 2) - hc011.*PD(2, 4).*PI(1, 2) - ...
   f011.*PD(2, 5).*PI(1, 2) - fc011.*PD(2, 6).*PI(1, 2) - h011.*PD(3, 3).*PI(1, 3) - hc011.*PD(3, 4).*PI(1, 3) - f011.*PD(3, 5).*PI(1, 3) - fc011.*PD(3, 6).*PI(1, 3) - ...
   (3.*gamma.*xi1.^2.*P(1, 2).*PI(1, 4))./m1 - (3.*gamma.*xi1.^2.*P(1, 2).*PI(1, 4))./Mf - (6.*gamma.*xi0.*xi2.*P(1, 2).*PI(1, 4))./m1 - (6.*gamma.*xi0.*xi2.*P(1, 2).*PI(1, 4))./Mf - ...
   (6.*gamma.*h011.*xi0.*xi1.*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h011.*xi0.*xi1.*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*hc011.*xi0.*xi1.*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc011.*xi0.*xi1.*P(1, 4).*PI(1, 4))./Mf - (6.*f011.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./m1 - (6.*f011.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*fc011.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./m1 - (6.*fc011.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./Mf - h011.*PD(4, 3).*PI(1, 4) - hc011.*PD(4, 4).*PI(1, 4) - f011.*PD(4, 5).*PI(1, 4) - ...
   fc011.*PD(4, 6).*PI(1, 4) + (3.*gamma.*xi1.^2.*P(1, 2).*PI(1, 5))./m1 + (6.*gamma.*xi0.*xi2.*P(1, 2).*PI(1, 5))./m1 + (6.*gamma.*h011.*xi0.*xi1.*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc011.*xi0.*xi1.*P(1, 4).*PI(1, 5))./m1 + (6.*f011.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 5))./m1 + (6.*fc011.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 5))./m1 - h011.*PD(5, 3).*PI(1, 5) - ...
   hc011.*PD(5, 4).*PI(1, 5) - f011.*PD(5, 5).*PI(1, 5) - fc011.*PD(5, 6).*PI(1, 5) - h011.*PD(6, 3).*PI(1, 6) - hc011.*PD(6, 4).*PI(1, 6) - f011.*PD(6, 5).*PI(1, 6) - ...
   fc011.*PD(6, 6).*PI(1, 6)) + epsilon.*u2.^2.*(-(h020.*PD(1, 3).*PI(1, 1)) - hc020.*PD(1, 4).*PI(1, 1) - f020.*PD(1, 5).*PI(1, 1) - fc020.*PD(1, 6).*PI(1, 1) - h020.*PD(2, 3).*PI(1, 2) - ...
   hc020.*PD(2, 4).*PI(1, 2) - f020.*PD(2, 5).*PI(1, 2) - fc020.*PD(2, 6).*PI(1, 2) - h020.*PD(3, 3).*PI(1, 3) - hc020.*PD(3, 4).*PI(1, 3) - f020.*PD(3, 5).*PI(1, 3) - ...
   fc020.*PD(3, 6).*PI(1, 3) - (3.*gamma.*xi1.*P(1, 2).^2.*PI(1, 4))./m1 - (3.*gamma.*xi1.*P(1, 2).^2.*PI(1, 4))./Mf - (6.*gamma.*h020.*xi0.*xi1.*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h020.*xi0.*xi1.*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h011.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h011.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc020.*xi0.*xi1.*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc020.*xi0.*xi1.*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc011.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc011.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f020.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./m1 - (6.*f020.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f011.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f011.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc020.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc020.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./Mf - (6.*fc011.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc011.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf - h020.*PD(4, 3).*PI(1, 4) - ...
   hc020.*PD(4, 4).*PI(1, 4) - f020.*PD(4, 5).*PI(1, 4) - fc020.*PD(4, 6).*PI(1, 4) + (3.*gamma.*xi1.*P(1, 2).^2.*PI(1, 5))./m1 + (6.*gamma.*h020.*xi0.*xi1.*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*h011.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*hc020.*xi0.*xi1.*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc011.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + ...
   (6.*f020.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 5))./m1 + (6.*f011.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc020.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 5))./m1 + ...
   (6.*fc011.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1 - h020.*PD(5, 3).*PI(1, 5) - hc020.*PD(5, 4).*PI(1, 5) - f020.*PD(5, 5).*PI(1, 5) - fc020.*PD(5, 6).*PI(1, 5) - ...
   h020.*PD(6, 3).*PI(1, 6) - hc020.*PD(6, 4).*PI(1, 6) - f020.*PD(6, 5).*PI(1, 6) - fc020.*PD(6, 6).*PI(1, 6)) + ...
 epsilon.^2.*u1.*(-(h101.*PD(1, 3).*PI(1, 1)) - hc101.*PD(1, 4).*PI(1, 1) - f101.*PD(1, 5).*PI(1, 1) - fc101.*PD(1, 6).*PI(1, 1) - h101.*PD(2, 3).*PI(1, 2) - hc101.*PD(2, 4).*PI(1, 2) - ...
   f101.*PD(2, 5).*PI(1, 2) - fc101.*PD(2, 6).*PI(1, 2) - h101.*PD(3, 3).*PI(1, 3) - hc101.*PD(3, 4).*PI(1, 3) - f101.*PD(3, 5).*PI(1, 3) - fc101.*PD(3, 6).*PI(1, 3) - ...
   (3.*gamma.*xi1.^2.*P(1, 1).*PI(1, 4))./m1 - (3.*gamma.*xi1.^2.*P(1, 1).*PI(1, 4))./Mf - (6.*gamma.*xi0.*xi2.*P(1, 1).*PI(1, 4))./m1 - (6.*gamma.*xi0.*xi2.*P(1, 1).*PI(1, 4))./Mf - ...
   (6.*gamma.*h101.*xi0.*xi1.*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h101.*xi0.*xi1.*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*hc101.*xi0.*xi1.*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc101.*xi0.*xi1.*P(1, 4).*PI(1, 4))./Mf - (6.*f101.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./m1 - (6.*f101.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*fc101.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./m1 - (6.*fc101.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./Mf - h101.*PD(4, 3).*PI(1, 4) - hc101.*PD(4, 4).*PI(1, 4) - f101.*PD(4, 5).*PI(1, 4) - ...
   fc101.*PD(4, 6).*PI(1, 4) + (3.*gamma.*xi1.^2.*P(1, 1).*PI(1, 5))./m1 + (6.*gamma.*xi0.*xi2.*P(1, 1).*PI(1, 5))./m1 + (6.*gamma.*h101.*xi0.*xi1.*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc101.*xi0.*xi1.*P(1, 4).*PI(1, 5))./m1 + (6.*f101.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 5))./m1 + (6.*fc101.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 5))./m1 - h101.*PD(5, 3).*PI(1, 5) - ...
   hc101.*PD(5, 4).*PI(1, 5) - f101.*PD(5, 5).*PI(1, 5) - fc101.*PD(5, 6).*PI(1, 5) - h101.*PD(6, 3).*PI(1, 6) - hc101.*PD(6, 4).*PI(1, 6) - f101.*PD(6, 5).*PI(1, 6) - ...
   fc101.*PD(6, 6).*PI(1, 6)) + epsilon.*u1.*u2.*(-(h110.*PD(1, 3).*PI(1, 1)) - hc110.*PD(1, 4).*PI(1, 1) - f110.*PD(1, 5).*PI(1, 1) - fc110.*PD(1, 6).*PI(1, 1) - h110.*PD(2, 3).*PI(1, 2) - ...
   hc110.*PD(2, 4).*PI(1, 2) - f110.*PD(2, 5).*PI(1, 2) - fc110.*PD(2, 6).*PI(1, 2) - h110.*PD(3, 3).*PI(1, 3) - hc110.*PD(3, 4).*PI(1, 3) - f110.*PD(3, 5).*PI(1, 3) - ...
   fc110.*PD(3, 6).*PI(1, 3) - (6.*gamma.*xi1.*P(1, 1).*P(1, 2).*PI(1, 4))./m1 - (6.*gamma.*xi1.*P(1, 1).*P(1, 2).*PI(1, 4))./Mf - (6.*gamma.*h110.*xi0.*xi1.*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h110.*xi0.*xi1.*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h011.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h011.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*h101.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h101.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*hc110.*xi0.*xi1.*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc110.*xi0.*xi1.*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc011.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc011.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc101.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc101.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f110.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./m1 - ...
   (6.*f110.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./Mf - (6.*f011.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f011.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f101.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f101.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc110.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc110.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./Mf - (6.*fc011.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - (6.*fc011.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - ...
   (6.*fc101.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc101.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf - h110.*PD(4, 3).*PI(1, 4) - hc110.*PD(4, 4).*PI(1, 4) - f110.*PD(4, 5).*PI(1, 4) - ...
   fc110.*PD(4, 6).*PI(1, 4) + (6.*gamma.*xi1.*P(1, 1).*P(1, 2).*PI(1, 5))./m1 + (6.*gamma.*h110.*xi0.*xi1.*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h011.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*h101.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*hc110.*xi0.*xi1.*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc011.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc101.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + (6.*f110.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 5))./m1 + (6.*f011.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f101.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc110.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 5))./m1 + (6.*fc011.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 + ...
   (6.*fc101.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1 - h110.*PD(5, 3).*PI(1, 5) - hc110.*PD(5, 4).*PI(1, 5) - f110.*PD(5, 5).*PI(1, 5) - fc110.*PD(5, 6).*PI(1, 5) - ...
   h110.*PD(6, 3).*PI(1, 6) - hc110.*PD(6, 4).*PI(1, 6) - f110.*PD(6, 5).*PI(1, 6) - fc110.*PD(6, 6).*PI(1, 6)) + ...
 epsilon.*u1.^2.*(-(h200.*PD(1, 3).*PI(1, 1)) - hc200.*PD(1, 4).*PI(1, 1) - f200.*PD(1, 5).*PI(1, 1) - fc200.*PD(1, 6).*PI(1, 1) - h200.*PD(2, 3).*PI(1, 2) - hc200.*PD(2, 4).*PI(1, 2) - ...
   f200.*PD(2, 5).*PI(1, 2) - fc200.*PD(2, 6).*PI(1, 2) - h200.*PD(3, 3).*PI(1, 3) - hc200.*PD(3, 4).*PI(1, 3) - f200.*PD(3, 5).*PI(1, 3) - fc200.*PD(3, 6).*PI(1, 3) - ...
   (3.*gamma.*xi1.*P(1, 1).^2.*PI(1, 4))./m1 - (3.*gamma.*xi1.*P(1, 1).^2.*PI(1, 4))./Mf - (6.*gamma.*h200.*xi0.*xi1.*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h200.*xi0.*xi1.*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*h101.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h101.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*hc200.*xi0.*xi1.*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc200.*xi0.*xi1.*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc101.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc101.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - ...
   (6.*f200.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./m1 - (6.*f200.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 4))./Mf - (6.*f101.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - ...
   (6.*f101.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - (6.*fc200.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./m1 - (6.*fc200.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 4))./Mf - ...
   (6.*fc101.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - (6.*fc101.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - h200.*PD(4, 3).*PI(1, 4) - hc200.*PD(4, 4).*PI(1, 4) - f200.*PD(4, 5).*PI(1, 4) - ...
   fc200.*PD(4, 6).*PI(1, 4) + (3.*gamma.*xi1.*P(1, 1).^2.*PI(1, 5))./m1 + (6.*gamma.*h200.*xi0.*xi1.*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h101.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc200.*xi0.*xi1.*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc101.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*f200.*gamma.*xi0.*xi1.*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f101.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + (6.*fc200.*gamma.*xi0.*xi1.*P(1, 6).*PI(1, 5))./m1 + (6.*fc101.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 - h200.*PD(5, 3).*PI(1, 5) - ...
   hc200.*PD(5, 4).*PI(1, 5) - f200.*PD(5, 5).*PI(1, 5) - fc200.*PD(5, 6).*PI(1, 5) - h200.*PD(6, 3).*PI(1, 6) - hc200.*PD(6, 4).*PI(1, 6) - f200.*PD(6, 5).*PI(1, 6) - ...
   fc200.*PD(6, 6).*PI(1, 6));

elseif order == 1
   du =l1.*u1  + ...
 u1.^2.*((-3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 4))./m1 - (3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 4))./Mf + (3.*gamma.*xi0.*P(1, 1).^2.*PI(1, 5))./m1) + ...
 u1.*u2.*((-6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 4))./m1 - (6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 4))./Mf + (6.*gamma.*xi0.*P(1, 1).*P(1, 2).*PI(1, 5))./m1) + ...
 u2.^2.*((-3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 4))./m1 - (3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 4))./Mf + (3.*gamma.*xi0.*P(1, 2).^2.*PI(1, 5))./m1) + ...
 u2.^3.*(-((gamma.*P(1, 2).^3.*PI(1, 4))./m1) - (gamma.*P(1, 2).^3.*PI(1, 4))./Mf - (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - ...
   (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (gamma.*P(1, 2).^3.*PI(1, 5))./m1 + (6.*gamma.*h020.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*hc020.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + ...
   (6.*f020.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc020.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1) + ...
  u1.*u2.^2.*((-3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 4))./m1 - (3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 4))./Mf - (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (3.*gamma.*P(1, 1).*P(1, 2).^2.*PI(1, 5))./m1 + (6.*gamma.*h020.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h110.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc020.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc110.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + (6.*f020.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f110.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc020.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 + (6.*fc110.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1) + ...
 u1.^2.*u2.*((-3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 4))./m1 - (3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 4))./Mf - (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./m1 - (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 4))./Mf - ...
   (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./m1 - ...
   (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 4))./Mf - (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - ...
   (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./m1 - (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 4))./Mf - (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf - (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./m1 - (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 4))./Mf + ...
   (3.*gamma.*P(1, 1).^2.*P(1, 2).*PI(1, 5))./m1 + (6.*gamma.*h110.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + (6.*gamma.*h200.*xi0.*P(1, 2).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc110.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*gamma.*hc200.*xi0.*P(1, 2).*P(1, 4).*PI(1, 5))./m1 + (6.*f110.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + ...
   (6.*f200.*gamma.*xi0.*P(1, 2).*P(1, 5).*PI(1, 5))./m1 + (6.*fc110.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1 + (6.*fc200.*gamma.*xi0.*P(1, 2).*P(1, 6).*PI(1, 5))./m1)+ u1.^3.*(-((gamma.*P(1, 1).^3.*PI(1, 4))./m1) - (gamma.*P(1, 1).^3.*PI(1, 4))./Mf - (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./m1 - ...
   (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 4))./Mf - (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./m1 - (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 4))./Mf - ...
   (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./m1 - (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 4))./Mf - (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./m1 - ...
   (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 4))./Mf + (gamma.*P(1, 1).^3.*PI(1, 5))./m1 + (6.*gamma.*h200.*xi0.*P(1, 1).*P(1, 3).*PI(1, 5))./m1 + ...
   (6.*gamma.*hc200.*xi0.*P(1, 1).*P(1, 4).*PI(1, 5))./m1 + (6.*f200.*gamma.*xi0.*P(1, 1).*P(1, 5).*PI(1, 5))./m1 + (6.*fc200.*gamma.*xi0.*P(1, 1).*P(1, 6).*PI(1, 5))./m1);

end    
Np = 1;
dy(1:Np,1) = real(du);
dy(Np+1:2*Np,1) = imag(du);

end