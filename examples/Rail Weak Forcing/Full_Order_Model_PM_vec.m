function dy = Full_Order_Model_PM_vec(t,y,F,epsilon,c,k,a,g,m,M,cf,Scal)

K = [a^2*(1 + m/M)*g*4/(Scal) - k*m/(M*(m + M)), k/M; k*m/(M + m)^2, -k/(M + m)];
Cd = [-cf*(M + m)/(M*m) - c*m/((M + m)*M), c/M; c*m/(m + M)^2, -c/(M + m)];

A = [zeros(2),eye(2);K,Cd];

Np = max(size(y))/4;
LinPart = A*[y(1:Np,1).';y(Np+1:2*Np,1).';y(2*Np+1:3*Np,1).';y(3*Np+1:4*Np,1).'];

dy(1:Np,1) = LinPart(1,:).';
dy(Np+1:2*Np,1) = LinPart(2,:).'; 
dy(2*Np+1:3*Np,1) = LinPart(3,:).' -4*g*(1+m/M)*y(1:Np,1).^3/Scal-16*(1+m/M)*a^4*y(1:Np,1).*y(2*Np+1:3*Np,1).^2/Scal^2;
dy(3*Np+1:4*Np,1) = LinPart(4,:).' +epsilon*F(t);
end