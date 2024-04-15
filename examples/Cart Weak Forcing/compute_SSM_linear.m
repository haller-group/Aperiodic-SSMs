function [y,modal] = compute_SSM_linear(xi_full_eval5,xi_full_eval3,xi_full_eval,t,lf,ls,gamma,V,VI,z,h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon)

% H2 =H_Coeff2(t);
% 
% H2 = H2.';
%     
% h401 = H2(1,:);
% h041 = H2(2,:);
% h311= H2(3,:);
% h131= H2(4,:);
% h221= H2(5,:);
% h302= H2(6,:);
% h032= H2(7,:);
% h212= H2(8,:);
% h122= H2(9,:);
% h203= H2(10,:);
% h023= H2(11,:);
% h113= H2(12,:);
% h104= H2(13,:);
% h014= H2(14,:);
% 
% 
% 
% F2 =F_Coeff2(t);
% F2 = F2.';
% 
% f401 = F2(1,:);
% f041 = F2(2,:);
% f311= F2(3,:);
% f131= F2(4,:);
% f221= F2(5,:);
% f302= F2(6,:);
% f032= F2(7,:);
% f212= F2(8,:);
% f122= F2(9,:);
% f203= F2(10,:);
% f023= F2(11,:);
% f113= F2(12,:);
% f104= F2(13,:);
% f014= F2(14,:);
% 
% 
% 
% H1 =H_Coeff1(t);
% 
% H1 = H1.';
% 
% F1 =F_Coeff1(t);
% 
% F1 = F1.';
% 
% h021 = H1(1,:);
% h201 = H1(2,:);
% h111 = H1(3,:);
% h012 = H1(4,:);
% h102 = H1(5,:);
% 
% f021 = F1(1,:);
% f201 = F1(2,:);
% f111 = F1(3,:);
% f012 = F1(4,:);
% f102 = F1(5,:);
% 
 yuk = conj(z);
% 
% p = h300*z.^3+h030*yuk.^3+h210*z.^2.*yuk+h120.*yuk.^2.*z ;
% p = p + h500*z.^5 + h050*yuk.^5 + h410*z.^4.*yuk + h140*z.*yuk.^4 + h230*z.^2.*yuk.^3 + h320*z.^3.*yuk.^2;
% p = p + epsilon^2*(h012.*yuk + h102.*z +h032.*yuk.^3 + h122.*z.*yuk.^2 + h212.*z.^2.*yuk + h302.*z.^3 );
% p = p + epsilon*(h021.*yuk.^2 + h201.*z.^2 + h111.*yuk.*z);
% 
% p = p + epsilon*(h041.*yuk.^4+h401.*z.^4 + h311.*z.^3.*yuk + h131.*z.*yuk.^3 + h221.*z.^2.*yuk.^2);
% p = p + epsilon^3*(h023.*yuk.^2 + h203.*z.^2 + h113.*z.*yuk);
% p = p + epsilon^4*(h014.*yuk + h104.*z);
% 
% r = f300*z.^3+f030*yuk.^3+f210*z.^2.*yuk+f120.*yuk.^2.*z ;
% r = r + f500*z.^5 + f050*yuk.^5 + f410*z.^4.*yuk + f140*z.*yuk.^4 + f230*z.^2.*yuk.^3 + f320*z.^3.*yuk.^2;
% r = r + epsilon^2*(f012.*yuk + f102.*z +f032.*yuk.^3 + f122.*z.*yuk.^2 + f212.*z.^2.*yuk + f302.*z.^3 );
% r = r + epsilon*(f021.*yuk.^2 + f201.*z.^2 + f111.*yuk.*z);
% 
% r = r + epsilon*(f041.*yuk.^4+f401.*z.^4 + f311.*z.^3.*yuk + f131.*z.*yuk.^3 + f221.*z.^2.*yuk.^2);
% r = r + epsilon^3*(f023.*yuk.^2 + f203.*z.^2 + f113.*z.*yuk);
% r = r + epsilon^4*(f014.*yuk + f104.*z);




vexc = [z;yuk;z.*0;z.*0;z.*0;z.*0]; % linear approximation
Solution = epsilon*xi_full_eval(t);%+epsilon^3*xi_full_eval3(t)+epsilon^5*xi_full_eval5(t);

TY = VI*Solution.';
% TU = [TY(3:end,:)*0;TY(3:end,:)];
y = real((V*vexc)) + Solution.';
modal = vexc + TY;
% y = real((V*vexc));


end