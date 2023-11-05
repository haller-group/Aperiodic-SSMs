hFig = figure('DefaultAxesFontSize',18);                       % Bring up new figure
% imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
% hFig.WindowState = 'maximized';
% pause(0.5)
delta = 0.5;
% tl =  100+2550;
tl = 50;
clear movieVector
shift_traj = y;
shift_trajZ = [real(modal(1,:));imag(modal(1,:));imag(modal(4,:));real(modal(3,:));real(modal(5,:))];%y;
ind = 1;
SP_TrajZ = inv(V)*SP_Traj.';
SP_TrajZ = [real(SP_TrajZ(1,:));imag(SP_TrajZ(1,:));imag(SP_TrajZ(4,:));real(SP_TrajZ(3,:));real(modal(5,:))];
SolutionZ = inv(V)*SolutionA.';
SolutionZ = [real(SolutionZ(1,:));imag(SolutionZ(1,:));imag(SolutionZ(4,:));real(SolutionZ(3,:));real(modal(5,:))];

Solution1 = SolutionA.';
SP_Traj1 = SP_Traj.';
nssm = 20;

for j = 150:2700 %2550:2550 %1000:1000 %150:150 %150:2700%2550 %1:20:5000%(max(size(t_sol))-tl)
i = j+0;

     clf
% subplot(2, 2, [1 3])
     rhosamp = linspace(0,0.5,201);
     plotdofs = [1 2 4]; 
     [RHO,Theta]=meshgrid(linspace(0,1,201),linspace(0,2*pi,201));
     x = linspace(-2,2,nssm);
     [X,Y]=meshgrid(x);
     Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*t_sol(j+tl);
     [ZZ,modalZ] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_ssm,l1,l1,gamma,V,inv(V),Z(:).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);
    
     %Z1 = reshape(real(modalZ(1,:)),201,201);
     %Z2 = reshape(imag(modalZ(1,:)),201,201);
     %Z3 = reshape(real(modalZ(4,:)),201,201);
     Z1 = reshape(ZZ(1,:),nssm,nssm);
     Z2 = reshape(ZZ(2,:),nssm,nssm);
     Z3 = reshape(ZZ(3,:),nssm,nssm);

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
     Z1 = reshape(real(modalZ(1,:)),nssm,nssm);
     Z2 = reshape(imag(modalZ(1,:)),nssm,nssm);
     Z3 = reshape(real(modalZ(3,:)),nssm,nssm);
%      Z1 = reshape(ZZ(1,:),201,201);
%      Z2 = reshape(ZZ(2,:),201,201);
%      Z3 = reshape(ZZ(3,:),201,201);

     h = surf(Z1,Z2,Z3);
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
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
  
%      plot3(SolutionZ(plotdofs(1),2000:2600),SolutionZ(plotdofs(2),2000:2600),SolutionZ(plotdofs(3),2000:2600),'.','LineWidth',3,'color',[0,0,1,0.5])
     hold on
     plot3(SP_TrajZ(plotdofs(1),j:i+tl),SP_TrajZ(plotdofs(2),j:i+tl),SP_TrajZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'Color',[0 0 0 0.5]);
     hold on
     plot3(shift_trajZ(plotdofs(1),j:i+tl),shift_trajZ(plotdofs(2),j:i+tl),shift_trajZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'Color',[1,0,0,0.5]);
     hold on
     plot3(SolutionZ(plotdofs(1),j:i+tl),SolutionZ(plotdofs(2),j:i+tl),SolutionZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'color',[0 0 1 0.5])
     hold on
%       plot3(SolutionZ(plotdofs(1),2000:j+tl),SolutionZ(plotdofs(2),2000:j+tl),SolutionZ(plotdofs(3),2000:j+tl),'-','LineWidth',3,'color',[0 1 1 0.5])
%      hold on
      plot3(shift_trajZ(plotdofs(1),i+tl),shift_trajZ(plotdofs(2),i+tl),shift_trajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color',[1,0,0]);
     hold on
     plot3(SP_TrajZ(plotdofs(1),i+tl),SP_TrajZ(plotdofs(2),i+tl),SP_TrajZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color','black');
     hold on
     plot3(SolutionZ(plotdofs(1),i+tl),SolutionZ(plotdofs(2),i+tl),SolutionZ(plotdofs(3),i+tl),'.','MarkerSize',30,'Color','blue');
     hold on
        
    xl = xlabel('Re[$u_1+\xi_{u_1}$]','Interpreter','latex');
     xl.Position = [-0.3395+2 -2.4427 -2.1373-0.3];
    xl.HorizontalAlignment = 'center';
    yl = ylabel('Im[$u_1+\xi_{u_1}$]','Interpreter','latex');
     yl.Position = [-2.3392 -0.6379+1.3 -1.9992-0.1];
    yl.HorizontalAlignment = 'center';
    zl = zlabel('Re[$v_2+\xi_{v_2}$]','Interpreter','latex');
    yticks([-2 0 2])

     
%      set(get(gca, 'ylabel' ), 'rotation' ,0,'VerticalAlignment' , 'top' , 'HorizontalAlignment' , 'center')
%      set(get(gca, 'xlabel' ), 'rotation' ,0)
%      set(get(gca, 'zlabel' ), 'rotation' ,0,'VerticalAlignment' , 'top','HorizontalAlignment' , 'center')
%      title('Modal Coordinates SSM for \epsilon = 0.1 and |\epsilon F_{max}| = 10')
%     legend('SSM - $W_{\epsilon}(x_{\epsilon}(t))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 5)$','Synchronised Anchor Trajectory region','Full Order Model','ROM - Prediction','$O(\epsilon^5)$ Anchor Trajectory ($x_{\epsilon}(t)$)','Interpreter','latex','location','northeast')
      legend('SSM - $W_{\epsilon}(x_{\epsilon}(t))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 5)$','Full order model','Reduced order model','Chaotic anchor trajectory ($x_{\epsilon}(t)$)','Interpreter','latex','location','northeast')

%     title('$\epsilon = 0.5$','Interpreter','latex');
    string = strcat('$t = $',num2str(round(t_sol(j+tl),2),'%4.2f' ),', max$(|F_{ext}|) = 3$ [$\mathbf{N}$]');
    if (j+tl) >=3000
     string = strcat('{\color{blue}t = ',num2str(round(t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.5');
     title(string);
    else
     title(string,'Interpreter','latex');
    end    
%     
%     axis([-0.5,0.5,-0.5,0.5,-0.02,0.02])
 view(-58,22)
    grid on
    box on 
     daspect([1,1,1])
   axis([-2.5,2.5,-2.5,2.5,-1.5,1.5])
%    axis([-3/2,3/2,-3/2,3/2,-3/2,3/2])
set(gcf,'color','white')

set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')

    figssm = gcf;

 movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

end

% myWriter = VideoWriter('Shake_Cart_SSM_video_fast','MPEG-4');
% myWriter.FrameRate = 50;
% 
% open(myWriter);
% 
% for indj = 1:2551
% writeVideo(myWriter,movieVector(indj));
% end
% 
% close(myWriter);

% set(hFig, 'Units' , 'Inches' );
% pos = get(hFig, 'Position' );
% set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(hFig, 'ssm_shaking_cart_snap3.pdf' , '-dpdf' , '-r300' )