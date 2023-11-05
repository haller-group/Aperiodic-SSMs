z = 0.5*exp(1i*0.6);
t = 0;
Omega = 0.01;
Opema =1.5; 
b = 2.2;
epk = 0;
IC_f = compute_SSM_phy(t,Omega,z,epk,Opema,b);
tf = 2*pi/(Omega)*4;
tspan = linspace(0,tf,5000);

ode_full = @(t,y) actual_model(t,y,Omega,b,epk,Opema);

tic
[tk,yy]=ode45(ode_full,tspan,real(IC_f));
toc


ode_reduced = @(t,zk) rom_temp_model(t,Omega,zk,epk,Opema,b);

tic
[tk,zz]=ode45(ode_reduced,tspan,[real(z);imag(z)]);
toc

% reduced to physical
klz = zz(:,1) + 1i*zz(:,2);
phy = [];
for i = 1:max(size(tspan))
    yl = compute_SSM_phy(tk(i),Omega,klz(i),epk,Opema,b);
    phy = [phy;real(yl.')];
end   

% Plot SSMs along trajectory - video
hFig = figure('DefaultAxesFontSize',20);                       % Bring up new figure
imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
hFig.WindowState = 'maximized';
pause(0.5)
delta = 0.5;
tl = 10;
clear movieVector
for i = 1:1%990%(max(size(traj.time))-tl)
 clf
 rhosamp = [0 0.01:0.01:2];
 [RHO,THETA] = meshgrid(rhosamp,0:2*pi/128:2*pi);
 zm = RHO.*exp(1i*THETA);
 SSMs = [];
  for ikl = 1:129
      for jkl = 1:201
        yml = compute_SSM_phy(tk(i+tl),Omega,zm(ikl,jkl),epk,Opema,b);
        SSMs = [SSMs;real(yml.')];   
      end
  end    
  h = surf(SSMs(:,1),SSMs(:,2),SSMs(:,3));

  h.EdgeColor = 'none';
  h.FaceColor = 'green';
  h.FaceAlpha = 0.3;

  view([1,1,1])
  grid on
  set(gca,'LineWidth',1.2);
  set(gca,'FontSize',30);
  xlabel(['$q_{\mathrm{',num2str(1),'}}$'],'interpreter','latex');
  ylabel(['$q_{\mathrm{',num2str(2),'}}$'],'interpreter','latex');
  zlabel(['$q_{\mathrm{',num2str(3),'}}$'],'interpreter','latex');


end
%NMTE Error

% Diff = sum(sqrt(sum((phy - yy).^2,2)),1)/((max(size(tk)))*max(sqrt(sum((phy).^2,2))));


