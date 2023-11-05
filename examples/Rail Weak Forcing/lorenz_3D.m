function dy = lorenz_3D(t,y)
Np = numel(y)/3;
sigma = 10;
rho = 28;
beta = 8/3;
dy(1:Np,1) =sigma*(y(Np+1:2*Np,1)-y(1:Np,1));
dy(Np+1:2*Np,1) = y(1:Np,1)*(rho - y(2*Np+1:3*Np,1))-y(Np+1:2*Np,1);
dy(2*Np+1:3*Np,1) = y(1:Np,1)*y(Np+1:2*Np,1)-beta*y(2*Np+1:3*Np,1);
end