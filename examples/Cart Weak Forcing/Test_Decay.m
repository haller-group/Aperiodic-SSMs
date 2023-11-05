ctspan = linspace(0,1000,8000);

m1 = 1;
m2 = 1;
Mf = 4;
kf = 6;
k = 1;
c = 0.3;
cf = 1.8;
gamma = 0.5;

FullS = @(t,y) full_system(t,y,epsilon,Force_Lorenz,m1,m2,Mf,kf,k,c,cf,gamma);

tic
IC = y0;
[tSP,SP_Traj] = ode45(FullS,ctspan,IC);
toc
y = y.';
figure 
indexR = 4;
hold on 
plot(ctspan,SP_Traj(:,indexR),'-','LineWidth',3,'color','red')
hold on 
