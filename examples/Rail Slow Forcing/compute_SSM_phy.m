function [Y,Modal,Net_Sol] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,t,Z,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha)
Y = [];
Modal = [];
Net_Sol = [];
for indf = 1:max(size(t))
alpha = epsilon*t(indf);
P = Valpha(Force_Lorenz(alpha),V);
PI = inv(P);
TempCell = num2cell(SSM_Coeff_A_1(alpha));
[h101,h200,h011,h110,h020,h300,h030,h120,h210]=deal(TempCell{:});
TempCell = num2cell(conj([h101,h200,h011,h110,h020,h300,h030,h120,h210]));
[hc101,hc200,hc011,hc110,hc020,hc300,hc030,hc120,hc210] = deal(TempCell{:});
TempCell = num2cell(SSM_Coeff_A_2(alpha));    
[h012,h102,h201,h021,h111] = deal(TempCell{:});

z = real(Z(indf));
yuk = imag(Z(indf));

p = h300*z.^3+h030*yuk.^3+h210*z.^2.*yuk+h120.*yuk.^2.*z ;
p = p + h200*z.^2 + h020*yuk.^2+h110*yuk*z;
p = p + epsilon^2*(h012.*yuk + h102.*z);
p = p + epsilon*(h021.*yuk.^2 + h201.*z.^2 + h111.*yuk.*z + h101*z + h011*yuk);



vexc = [z;yuk;p;conj(p)];

Solution = XI_0(alpha)+ epsilon*XI_1(alpha)+epsilon^2*XI_2(alpha)+epsilon^3*XI_3(alpha);

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