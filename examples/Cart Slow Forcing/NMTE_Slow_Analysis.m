
Mat_NMTE = zeros(3,2);
epsilonq = [0.001,0.008,0.01];
Orderq = [1,3];
zo = [0.3*exp(-1i*pi/2),0.3*exp(1i*pi/2),0.5*exp(-1i*pi/2),0.5*exp(1i*pi/2),0.4*exp(1i*0.5),0.4*exp(-1i*0.5)];
for indg = 1:2 
    order = Orderq(indg);
for indf = 1:3
    epsilon = epsilonq(indf);
    NMTET = [];
for ind = 1:6

ctspan = linspace(0,6,100000)/epsilon;
ROM=@(t,z) rom_temp_model_adiabatic_order(order,t,z,SSM_Coeff_A_2,SSM_Coeff_A_1,xi_01,xi_11,xi_21,Valpha,V,A,Force_Lorenz,Dalpha,gamma,epsilon);
q0 = zo(ind);
[y0,model0,Net_Sol0] = compute_SSM_phy_Order(order,XI_0,XI_1,XI_2,XI_3,ctspan(1),gamma,q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

tic 
[t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
toc

[y,modal,Net_Sol] = compute_SSM_phy_Order(order,XI_0,XI_1,XI_2,XI_3,t_sol,gamma,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz);
tic
IC =y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc


NMTE = sum(sqrt(sum((y.' - SP_Traj).^2,2)))/(sqrt(max(sum((SP_Traj).^2,2)))*max(size(ctspan)));
NMTET = [NMTET,NMTE];
end

% origUnits = fig.Units;
% fig.Units = fig.PaperUnits;
% fig.PaperSize = fig.Position(3:4);
% fig.Units = origUnits;
% exportgraphics(fig, 'e008xcdot_slow.pdf');

NMTEavg = sum(NMTET)/6;
Mat_NMTE(indf,indg) = NMTEavg;
end
end
save('NMTE_slow.mat',"Mat_NMTE")
