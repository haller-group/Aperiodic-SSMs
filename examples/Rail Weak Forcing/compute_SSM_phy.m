function [y,modal] = compute_SSM_phy(XI_hyperbolic_solution,t,V,VI,u,h030,h300,h120,h210,h111,h201,h021,h012,h102,epsilon)
u1 = u(1,:);
u2 = u(2,:);


v1 = h030.*u2.^3 + h300.*u1.^3 + h210.*u1.^2.*u2 + h120.*u1.*u2.^2;
% v1 = v1 + h500*u1.^5 + h050*u2.^5 + h410*u1.^4.*u2 + h140*u1.*u2.^4 + h230*u1.^2.*u2.^3 + h320*u1.^3.*u2.^2;

v1 = v1 + epsilon*(h021(t).*u2.^2 + h201(t).*u1.^2 + h111(t).*u1.*u2);
v1 = v1 + epsilon^2*(h012(t).*u2 +h102(t).*u1);

vexc = [u1;u2;v1;conj(v1)];
Solution = epsilon*XI_hyperbolic_solution{1,1}(t)+epsilon^3*XI_hyperbolic_solution{1,3}(t);%+epsilon^9*XI_hyperbolic_solution{1,9}(t);

TY = VI*Solution.';
% TU = [TY(3:end,:)*0;TY(3:end,:)];
% max(size(vexc(:,1)))
% kj = real((V*vexc))
y = real((V*vexc)) + Solution.';
modal = vexc + TY;
% y = real((V*vexc));


end