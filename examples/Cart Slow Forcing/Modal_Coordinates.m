hFig = figure('DefaultAxesFontSize',18);                       % Bring up new figure
% imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
% hFig.WindowState = 'maximized';
% pause(0.5)
delta = 0.5;
tl =  200-1;
% tl = 1000;
clear movieVector
shift_traj = y;
% shift_trajZ = [real(modal(1,:));imag(modal(1,:));imag(modal(4,:));real(modal(4,:))];%y;
% ind = 1;
% SP_TrajZ = inv(V)*SP_Traj.';
% SP_TrajZ = [real(SP_TrajZ(1,:));imag(SP_TrajZ(1,:));imag(SP_TrajZ(4,:));real(SP_TrajZ(4,:))];
% SolutionZ = inv(V)*SolutionA.';
% SolutionZ = [real(SolutionZ(1,:));imag(SolutionZ(1,:));imag(SolutionZ(4,:));real(SolutionZ(4,:))];
% 
% Solution1 = SolutionA.';
SP_Traj1 = SP_Traj.';

ind = 1;
SolutionA = Net_Sol;
nssm = 20;
for j = 1:1:2500 %1:1 %3500:3500 %500:500 %1:20:10000-tl % 1:100:70000%:1200%(max(size(t_sol))-tl)
i = 1+0;
     clf
% subplot(2, 2, [1 3])
     rhosamp = linspace(0,0.5,201);
     plotdofs = [3 6 2]; 
     [RHO,Theta]=meshgrid(linspace(0,2,nssm),linspace(0,2*pi,nssm));
     Z = RHO.*exp(1i*Theta);
     x = linspace(-1,1,nssm);
     [X,Y]=meshgrid(x);
     Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*tSP(j+tl);
%      [ZZ,modalZ] = compute_SSM_phy(xi_full_eval3,xi_full_eval,t_ssm,lf,ls,gamma,V,inv(V),Z(:).',h030,h300,h210,h120,h111,h021,h201,h102,h012,f030,f300,f210,f120,f111,f021,f201,f102,f012,epsilon);
[ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0,XI_1,XI_2,XI_3,t_ssm,gamma,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

     %Z1 = reshape(real(modalZ(1,:)),201,201);
     %Z2 = reshape(imag(modalZ(1,:)),201,201);
     %Z3 = reshape(real(modalZ(4,:)),201,201);
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);

%      h = surf(Z1,Z2,Z3)
%      h.EdgeColor = 'none';
%      h.FaceColor = [.7 .7 .7];
%      h.FaceAlpha = 0.3;
%      hold on 

%      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
     starto = 1;
%      plot3(Solution1(plotdofs(1),starto:end),Solution1(plotdofs(2),starto:end),Solution1(plotdofs(3),starto:end),'--','LineWidth',1,'color','green')
%      hold on 
%      plot3(Solution1(plotdofs(1),starto),Solution1(plotdofs(2),starto),Solution1(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(Solution1(plotdofs(1),end),Solution1(plotdofs(2),end),Solution1(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
% 
%      plot3(SP_Traj(i:i+tl,plotdofs(1)),SP_Traj(i:i+tl,plotdofs(2)),SP_Traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color','black');
%      hold on
%      plot3(shift_traj(plotdofs(1),i:i+tl),shift_traj(plotdofs(2),i:i+tl),shift_traj(plotdofs(3),i:i+tl),'--','LineWidth',3,'Color','red');
%      hold on
%      plot3(Solution1(plotdofs(1),i:i+tl),Solution1(plotdofs(2),i:i+tl),Solution1(plotdofs(3),i:i+tl),'-o','LineWidth',3,'color','green')
%      hold on
% 
%      plot3(SP_Traj(i+tl,plotdofs(1)),SP_Traj(i+tl,plotdofs(2)),SP_Traj(i+tl,plotdofs(3)),'.','MarkerSize',30,'Color','black');
%      hold on 
%      plot3(shift_traj(plotdofs(1),i+tl),shift_traj(plotdofs(2),i+tl),shift_traj(plotdofs(3),i+tl),'.','MarkerSize',30,'Color','red');%[0.4940 0.1840 0.5560]
%      hold on
%      plot3(Solution1(plotdofs(1),i+tl),Solution1(plotdofs(2),i+tl),Solution1(plotdofs(3),i+tl),'.','MarkerSize',30,'color','green')
%      hold on
% 
% %      plot3(shift_trajZ(plotdofs(1),i:i+tl),shift_trajZ(plotdofs(2),i:i+tl),shift_trajZ(plotdofs(3),i:i+tl),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %      hold on
% %      plot3(shift_trajZ(plotdofs(1),i+tl),shift_trajZ(plotdofs(2),i+tl),shift_trajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %      hold on 
% %      plot3(SP_TrajZ(plotdofs(1),i:i+tl),SP_TrajZ(plotdofs(2),i:i+tl),SP_TrajZ(plotdofs(3),i:i+tl),'-','LineWidth',3,'Color','black');
% %      hold on
% %      plot3(SP_TrajZ(plotdofs(1),i+tl),SP_TrajZ(plotdofs(2),i+tl),SP_TrajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color','black');
% %     
%     
%      
%     legend('FOM','ROM','O(\epsilon) Anchor','Location' , 'best' )%,'FOM','FOM head')
%  
% %     h
% %     hold on
% 
% 
%      xlabel('$q_1 [m]$','Interpreter','latex');
%      ylabel('$q_2 [m]$','Interpreter','latex');
%      zlabel('$\dot{q}_2 [m/s]$','Interpreter','latex');
%      %   daspect([1,1,1])
% axis([-1,1,-1,1,-3,3])
%     view(32,26)
%     grid on
% %      W0 = WTF(target(Omega*traj.time(i+tl)));
% %   
% % coord_change_vec2 = 0;
% % coord_change_vec = [traj_actual(i+tl,plotdofs(1));traj_actual(i+tl,plotdofs(2));traj_actual(i+tl,plotdofs(3))];
% %      [zdof1,zdof2,zdof3]=plot_2D_auto_SSM(coord_change_vec2,coord_change_vec,W0,rhosamp,plotdofs,{'$q_1$ m','$q_2$ m','$\dot{q}_2$ m/s'});
% %      hold on 
% %      plot3(traj_actual1(:,plotdofs(1)),traj_actual1(:,plotdofs(2)),traj_actual1(:,plotdofs(3)),'-','Linewidth',3,'color','red');
% %     hold on
% %     plot3(coord_change_vec(1),coord_change_vec(2),coord_change_vec(3),'.','MarkerSize',30,'color','blue');
% %     hold on
% %     plot3(shift_traj(i:i+tl,plotdofs(1)),shift_traj(i:i+tl,plotdofs(2)),shift_traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(shift_traj(i+tl,plotdofs(1)),shift_traj(i+tl,plotdofs(2)),shift_traj(i+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(shift_traj(i:i+tl,plotdofs(1)),shift_traj(i:i+tl,plotdofs(2)),shift_traj(i:i+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0.4940 0.1840 0.5560]);
% %     hold on 
% %     plot3(SP_Traj_Off(i:i+2*tl,plotdofs(1)),SP_Traj_Off(i:i+2*tl,plotdofs(2)),SP_Traj_Off(i:i+2*tl,plotdofs(3)),'-','LineWidth',3,'Color','black');
% %     hold on 
% %     plot3(SP_Traj_Off(i+2*tl,plotdofs(1)),SP_Traj_Off(i+2*tl,plotdofs(2)),SP_Traj_Off(i+2*tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0.4940 0.1840 0.5560]);
% %     
% %     hold on 
% %     view([108 8])
% %     axis([-0.02,0.02,-0.02,0.02,-0.02,0.02])
% %    daspect([1 1 1])
% %     axis([-0.1,0.1,-0.1,0.1,-0.4,0.4])
%     
% title('Tracking ROM prediction, \epsilon = 0.5,  |\epsilon F|_{max} = 50')
% subplot(2, 2, 2)
% % rhosamp = linspace(0,2,201);
% %      plotdofs = [1 2 4]; 
% %      [RHO,Theta]=meshgrid(linspace(0,2,201),linspace(0,2*pi,201));
% %      Z = epsilon^(1/Order)*RHO.*exp(1i*Theta);
% %      t_ssm = ones(1,201*201)*t_sol(i+tl);
% %      [ZZ,modalZ] = compute_SSM_phy(xi_net_sol,t_ssm,lf,ls,gamma,V,inv(V),Z(:).');
%      %Z1 = reshape(real(modalZ(1,:)),201,201);
%      %Z2 = reshape(imag(modalZ(1,:)),201,201);
%      %Z3 = reshape(real(modalZ(4,:)),201,201);
%      Z1 = reshape(ZZ(1,:),201,201);
%      Z2 = reshape(ZZ(2,:),201,201);
%      Z3 = reshape(ZZ(3,:),201,201);
% 
%      h = surf(Z1,Z2,Z3)
%      h.EdgeColor = 'none';
%      h.FaceColor = [.7 .7 .7];
%      h.FaceAlpha = 0.3;
%      hold on 
% 
% %      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
% %      hold on 
% %      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
% %      hold on 
%      starto = 1;
%      plot3(Solution1(plotdofs(1),starto:end),Solution1(plotdofs(2),starto:end),Solution1(plotdofs(3),starto:end),'-o','LineWidth',3,'color','green')
%      hold on 
%      plot3(Solution1(plotdofs(1),starto),Solution1(plotdofs(2),starto),Solution1(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(Solution1(plotdofs(1),end),Solution1(plotdofs(2),end),Solution1(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
%      
%       plot3(SP_Traj(j:j+tl,plotdofs(1)),SP_Traj(j:j+tl,plotdofs(2)),SP_Traj(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color','black');
%      hold on
%      plot3(shift_traj(plotdofs(1),j:j+tl),shift_traj(plotdofs(2),j:j+tl),shift_traj(plotdofs(3),j:j+tl),'--','LineWidth',3,'Color','red');
%      hold on
%     
%      plot3(shift_traj(plotdofs(1),j+tl),shift_traj(plotdofs(2),j+tl),shift_traj(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','red');
%      hold on
%      plot3(SP_Traj(j+tl,plotdofs(1)),SP_Traj(j+tl,plotdofs(2)),SP_Traj(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color','black');
%      hold on
%      xlabel('$q_1 [m]$','Interpreter','latex');
%      ylabel('$q_2 [m]$','Interpreter','latex');
%      zlabel('$\dot{q}_2 [m/s]$','Interpreter','latex');
%        % daspect([1,1,1])
%      legend('SSM - W_{\epsilon}(x_{\epsilon}(t)) O(\epsilon^2 u^3, u^5,\epsilon)','Location' , 'best' )%,'FOM','FOM head')
%  view(-71,19)
%     grid on
% axis([-2,2,-2,2,-2,2])
% title('Physical Coordinates SSM')
% subplot(2, 2, 4)
% rhosamp = linspace(0,2,201);
%      plotdofs = [1 2 4]; 
%      [RHO,Theta]=meshgrid(linspace(0,2,201),linspace(0,2*pi,201));
%      Z = epsilon^(1/Order)*RHO.*exp(1i*Theta);
%      t_ssm = ones(1,201*201)*t_sol(i+tl);
%      [ZZ,modalZ] = compute_SSM_phy(xi_net_sol,t_ssm,lf,ls,gamma,V,inv(V),Z(:).');
%      Z1 = reshape(real(modalZ(1,:)),201,201);
%      Z2 = reshape(imag(modalZ(1,:)),201,201);
%      Z3 = reshape(real(modalZ(4,:)),201,201);
%      Z1 = reshape(ZZ(1,:),201,201);
%      Z2 = reshape(ZZ(2,:),201,201);
%      Z3 = reshape(ZZ(3,:),201,201);

     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 

%      plot3(SolutionZ(plotdofs(1),:),SolutionZ(plotdofs(2),:),SolutionZ(plotdofs(3),:),'-','LineWidth',3,'color','red')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),1),SolutionZ(plotdofs(2),1),SolutionZ(plotdofs(3),1),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%      hold on 
     starto = 1;
     
%      plot3(SolutionZ(plotdofs(1),starto),SolutionZ(plotdofs(2),starto),SolutionZ(plotdofs(3),starto),'.','MarkerSize',30,'color','blue')
%      hold on 
%      plot3(SolutionZ(plotdofs(1),end),SolutionZ(plotdofs(2),end),SolutionZ(plotdofs(3),end),'.','MarkerSize',30,'color','magenta')
%     hold on
%     plot3(SP_TrajZ(plotdofs(1),j:j+tl),SP_TrajZ(plotdofs(2),j:j+tl),SP_TrajZ(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color','black');
%      hold on

     plot3(SolutionA(plotdofs(1),1:1:end),SolutionA(plotdofs(2),1:1:end),SolutionA(plotdofs(3),1:1:end),'-.','color',[0,0,1,0.8])
     hold on
     flm = plot3(SP_Traj1(plotdofs(1),j:j+tl),SP_Traj1(plotdofs(2),j:j+tl),SP_Traj1(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color',[0 0 0 0.5]);
     hold on
     roml = plot3(shift_traj(plotdofs(1),j:j+tl),shift_traj(plotdofs(2),j:j+tl),shift_traj(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color',[1 0 0 0.5]);
     hold on
     anl = plot3(SolutionA(plotdofs(1),j:j+tl),SolutionA(plotdofs(2),j:j+tl),SolutionA(plotdofs(3),j:j+tl),'-','LineWidth',3,'color',[0,0,1,0.8]);
     hold on
     plot3(shift_traj(plotdofs(1),j+tl),shift_traj(plotdofs(2),j+tl),shift_traj(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','red');
     hold on
     plot3(SP_Traj1(plotdofs(1),j+tl),SP_Traj1(plotdofs(2),j+tl),SP_Traj1(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','black');
     hold on
     plot3(SolutionA(plotdofs(1),j+tl),SolutionA(plotdofs(2),j+tl),SolutionA(plotdofs(3),j+tl),'.','MarkerSize',30,'Color',[0,0,1]);
     hold on
      
    xl = xlabel('$x_c \,[$m$]$','Interpreter','latex');
    xl.Position = [0.2830+4 -0.6436 -1.3701-0.2];
    xl.HorizontalAlignment = 'center';
    xticks([-10 0 10])
    yl = ylabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    yl.Position = [-10.0940-0.2 -0.0183 -1.2803+0.3];
    yticks([-1 0 1])
     zlabel('$q_2 \,[$m$]$','Interpreter','latex');
     zticks([-0.5 0 0.5])
%        title('Adiabatic SSM reduction, $\epsilon = 0.008$','Interpreter','latex')
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
     lk =  legend([h flm roml anl],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Chaotic anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
 lk.Position = [0.0023 0.6272 0.5865 0.1986];
    string = strcat('$\alpha = $',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', $\epsilon = 0.008$');
    if (round(epsilon*t_sol(j+tl),2)) >=8
     string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
     title(string,'Interpreter','latex');
    else
     title(string,'Interpreter','latex');
    end   
   axis([-10,11,-1.2,1.2,-0.5,0.6])
   view(-79,7)
    grid on
    box on
           daspect([1,1,1])
%  axis([-0.2,0.2,-0.2,0.2,-0.01,0.01])
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

end

% myWriter = VideoWriter('SSM_Shaking_Cart_adiabatic_faster','MPEG-4');
% % myWriter.FrameRate = 40;
% % 
% % open(myWriter);
% % writeVideo(myWriter,movieVector);
% % close(myWriter);
% 
% myWriter.FrameRate = 80;
% 
% open(myWriter);
% 
% for indj = 1:2500
% writeVideo(myWriter,movieVector(indj));
% end
% 
% close(myWriter);

% set(hFig, 'Units' , 'Inches' );
% pos = get(hFig, 'Position' );
% set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(hFig, 'ssm_shaking_cart_slow_snap2.pdf' , '-dpdf')