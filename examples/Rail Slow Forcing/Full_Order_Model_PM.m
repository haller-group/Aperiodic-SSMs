function dy = Full_Order_Model_PM(t,y,F,epsilon,c,k,a,g,m,Mf,cf,scale)

K = -[-k,k*m/(m+Mf);k*m/(m+Mf),4*a^2*m*g/scale-k*m^2/(m+Mf)^2];
C = -[-c,c*m/(m+Mf);c*m/(m+Mf),-cf-c*m^2/(m+Mf)^2];
M = [m+Mf,0;0,m*Mf/(m+Mf)];

A = [zeros(2),eye(2);-inv(M)*K,-inv(M)*C];

LinPart = A*y;
Np =1;
dy(1:Np,1) = LinPart(1);
dy(Np+1:2*Np,1) = LinPart(2); 
dy(2*Np+1:3*Np,1) = LinPart(3) +F(epsilon*t)/(m+Mf);
dy(3*Np+1:4*Np,1) = LinPart(4) -4*g/scale*(1+m/Mf)*y(2)^3 -16*(1+m/Mf)*a^4/scale^2*y(2)*y(4)^2;
end