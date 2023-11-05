function [Y,Modal,Net_Sol] = compute_SSM_phy_Order(order,XI_0,XI_1,XI_2,XI_3,t,gamma,Z,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha)
Y = [];
Modal = [];
Net_Sol = [];
for indf = 1:max(size(t))
alpha = epsilon*t(indf);
P = Valpha(Force_Lorenz(alpha),V);
PI = inv(P);
TempCell = num2cell(SSM_Coeff_A_1(alpha));
[h101,h200,h011,h110,h020,f101,f200,f011,f110,f020,h300,h030,h120,h210,f300,f030,f120,f210]=deal(TempCell{:});
TempCell = num2cell(conj([h101,h200,h011,h110,h020,f101,f200,f011,f110,f020,h300,h030,h120,h210,f300,f030,f120,f210]));
[hc101,hc200,hc011,hc110,hc020,fc101,fc200,fc011,fc110,fc020,hc300,hc030,hc120,hc210,fc300,fc030,fc120,fc210] = deal(TempCell{:});
TempCell = num2cell(SSM_Coeff_A_2(alpha));    
[h012,h102,h201,h021,h111,f012,f102,f201,f021,f111] = deal(TempCell{:});
if order == 3
z = Z(indf);
yuk = conj(z);

p = h300*z.^3+h030*yuk.^3+h210*z.^2.*yuk+h120.*yuk.^2.*z ;
p = p + h200*z.^2 + h020*yuk.^2+h110*yuk*z;
p = p + epsilon^2*(h012.*yuk + h102.*z);
p = p + epsilon*(h021.*yuk.^2 + h201.*z.^2 + h111.*yuk.*z + h101*z + h011*yuk);

r = f300*z.^3+f030*yuk.^3+f210*z.^2.*yuk+f120.*yuk.^2.*z ;
r = r + f200*z.^2 + f020*yuk.^2+f110*yuk*z;
r = r + epsilon^2*(f012.*yuk + f102.*z);
r = r + epsilon*(f021.*yuk.^2 + f201.*z.^2 + f111.*yuk.*z + f101*z + f011*yuk);


vexc = [z;yuk;p;conj(p);r;conj(r)];

Solution = XI_0(alpha)+ epsilon*XI_1(alpha)+epsilon^2*XI_2(alpha)+epsilon^3*XI_3(alpha);

TY = PI*Solution.';
% TU = [TY(3:end,:)*0;TY(3:end,:)];
y = real((P*vexc)) + Solution.';
Y = [Y,y];
modal = vexc + TY;
Modal = [Modal,modal];
Net_Sol = [Net_Sol,Solution.'];
% y = real((V*vexc));
elseif order == 1
    z = Z(indf);
yuk = conj(z);

p = h300*z.^3+h030*yuk.^3+h210*z.^2.*yuk+h120.*yuk.^2.*z ;
p = p + h200*z.^2 + h020*yuk.^2+h110*yuk*z;

r = f300*z.^3+f030*yuk.^3+f210*z.^2.*yuk+f120.*yuk.^2.*z ;
r = r + f200*z.^2 + f020*yuk.^2+f110*yuk*z;


vexc = [z;yuk;p;conj(p);r;conj(r)];

Solution = XI_0(alpha)+ epsilon*XI_1(alpha);

TY = PI*Solution.';
% TU = [TY(3:end,:)*0;TY(3:end,:)];
y = real((P*vexc)) + Solution.';
Y = [Y,y];
modal = vexc + TY;
Modal = [Modal,modal];
Net_Sol = [Net_Sol,Solution.'];
% y = real((V*vexc));
end

end

end