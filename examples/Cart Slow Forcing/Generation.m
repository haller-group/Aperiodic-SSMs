% Plotting errors for variety of initial conditions 
rho = linspace(0.05,0.4,10);
theta = linspace(0,2*pi,10);
zq = rho.*exp(1i*theta);
combined_cell = cell(max(size(zq)),5);
figure
epsilon = 0.001;
for ij = 1:max(size(zq))

q0 = zq(ij);
ctspan = linspace(6,20,90000)/epsilon;
[y0,model0,Net_Sol0] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,ctspan(1),gamma,q0,epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

IC = [real(q0);imag(q0)];

tic 
[t_sol,ROM_sol] = ode45(ROM,ctspan,IC);
toc

[y,modal,Net_Sol] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,t_sol,gamma,(ROM_sol(:,1)+1i*ROM_sol(:,2)).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);


FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz);
tic
IC =y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc

% y = y.';
% tic
%  [F, lambda, V1, G, DG] = functionFromTensors(M, C, K, fnl);
% toc
% off = 0;
% Gh = @(t,x) F(t,x)+epsilon*[0;0;Force_Lorenz(t);Force_Lorenz(t)];
% Solution = epsilon*XI_hyperbolic_solution{1,1}(t_sol.');
% tic
% IC = y0;
% [tSP,SP_Traj] = ode45(Gh,ctspan,IC);
% toc


rel_error1 = sqrt(sum((y - Net_Sol).^2))/max(sqrt(sum((Net_Sol).^2)));

rel_error2 = sqrt(sum((y - SP_Traj.').^2))/max(sqrt(sum((SP_Traj.').^2)));
combined_cell{ij,1} = y;
combined_cell{ij,2} = SP_Traj.';
combined_cell{ij,3} = Net_Sol.';
combined_cell{ij,4} = rel_error1;
combined_cell{ij,5} = rel_error2;

y = y.';
SolutionA = Net_Sol.';
subplot(3,2,1)
indexR = 1;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$q_1 \,[$m$]$','Interpreter','latex');
legend('FOM','ROM $O(\epsilon^{k_1} u^{k_2}, k_1 + k_2 = 3)$','$O(\epsilon^3)$ Anchor Trajectory')
hold on 
subplot(3,2,2)
indexR = 2;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$q_2 \,[$m$]$','Interpreter','latex');
hold on

subplot(3,2,3)
indexR = 3;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$X_{CM} \,[$m$]$','Interpreter','latex');
hold on

subplot(3,2,4)
indexR = 4;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$\dot{q}_1 \,[$m/s$]$','Interpreter','latex');
hold on

subplot(3,2,5)
indexR = 5;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$\dot{q}_2 \,[$m/s$]$','Interpreter','latex');
hold on

subplot(3,2,6)
indexR = 6;
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','black')
hold on 
plot(ctspan,y(:,indexR),'--','LineWidth',3,'color','red')
hold on 
plot(ctspan,SolutionA(:,indexR),'--','LineWidth',3,'color',[0 1 0 0.3])

xlabel('$t \,[$s$]$','Interpreter','latex');
ylabel('$\dot{X}_{CM} \,[$m/s$]$','Interpreter','latex');
hold on

end

figure
l1 =0;
l2 = 0;
for ij =  1:max(size(zq))
    subplot(1,2,1)
    plot(ctspan,combined_cell{ij,4},'-','LineWidth',2,'color',[.7 .7 .7])
    hold on 
    xlabel('$t \,[$s$]$','Interpreter','latex');
    ylabel('Relative Error','Interpreter','latex');
    title('Relative Error between $O(\epsilon^9)$ target solution and ROM prediction')
    subplot(1,2,2)
    plot(ctspan,combined_cell{ij,5},'-','LineWidth',2,'color',[.7 .7 .7])
    hold on 
    xlabel('$t \,[$s$]$','Interpreter','latex');
    ylabel('Relative Error','Interpreter','latex');
    title('Relative Error between true solution and ROM prediction')
    l1 = l1 + combined_cell{ij,4};
    l2 = l2 + combined_cell{ij,5};
end    


subplot(1,2,1)
    plot(ctspan,l1/10,'-','LineWidth',2,'color','black')
    hold on 
    xlabel('$t \,[$s$]$','Interpreter','latex');
    ylabel('Relative Error','Interpreter','latex');
    title('Relative Error between $O(\epsilon)$ target solution and ROM prediction')
    subplot(1,2,2)
    plot(ctspan,l2/10,'-','LineWidth',2,'color','black')
    hold on 
    xlabel('$t \,[$s$]$','Interpreter','latex');
    ylabel('Relative Error','Interpreter','latex');
    title('Relative Error between true solution and ROM prediction')




