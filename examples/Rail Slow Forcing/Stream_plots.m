hFig = figure('DefaultAxesFontSize',20);                       % Bring up new figure
% imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
% hFig.WindowState = 'maximized';
pause(0.5)
delta = 0.5;
tl =  30;
% tl = 1000
clear movieVector

nssm = 20;
m0=1;
m1 = 1;

m2 = 1;
m3 = 1;
spand = 50;
ind = 1;
for marker = 2:4 %to(4)+10:to(4)+10%2:4 %to(4)+10:to(4)+10%1:spand:17000 % 1:100:70000%:1200%(max(size(t_sol))-tl)

     
% subplot(2, 2, [1 3])
   
     

       
%       plot3(store_Traj{1,3}(1:50:end,plotdofs(1)),store_Traj{1,3}(1:50:end,plotdofs(2)),store_Traj{1,3}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(1:50:end,plotdofs(1)),store_Traj{1,4}(1:50:end,plotdofs(2)),store_Traj{1,4}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
     for j = 1:40
         clf
         plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-','LineWidth',2,'Color',[1 0 0 0.5]);
       hold on 
         plotdofs = [2 4 3]; 
     m0 = j;
     x = linspace(-1,1,nssm);
     y = linspace(-0.5,0.5,nssm);
     [X,Y]=meshgrid(x,y);
     
    if m0+tl <=40
                  Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
          unst = plot3(store_Traj{1,2}(j+to(marker-1),plotdofs(1)),store_Traj{1,2}(j+to(marker-1),plotdofs(2)),store_Traj{1,2}(j+to(marker-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to(marker-1),plotdofs(1)),store_Traj{1,2}(j+to(marker-1),plotdofs(2)),store_Traj{1,2}(j+to(marker-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  

%               unst = plot3(store_Traj{1,2}(j+to(marker-1)-1:j+to(marker-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to(marker-1)-1:j+to(marker-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to(marker-1)-1:j+to(marker-1)-1+tl,plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
% %       hold on 
% %       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       
% %       hold on
%        plot3(store_Traj{1,2}(j+to(marker-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to(marker-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to(marker-1)-1+tl,plotdofs(3)),'.','MarkerSize',30,'Color','red');
%       hold on 
% %       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       
      elseif m0<=40
          end1=40;
                 Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
     unst = plot3(store_Traj{1,2}(j+to(marker-1),plotdofs(1)),store_Traj{1,2}(j+to(marker-1),plotdofs(2)),store_Traj{1,2}(j+to(marker-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to(marker-1),plotdofs(1)),store_Traj{1,2}(j+to(marker-1),plotdofs(2)),store_Traj{1,2}(j+to(marker-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  
    end

       for nTraj = 1:16
             
             if m0+tl <=40
                 
                 ft = plot3(store_Traj{marker,4}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
            rt = plot3((store_Traj{marker,2}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
             plot3((store_Traj{marker,2}(m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(3)+(nTraj-1)*4)),'.','MarkerSize',30,'Color','green');
            hold on
             plot3(store_Traj{marker,4}(m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','black');
            hold on
            elseif m0<=40
                
                end1=40;
                ft = plot3(store_Traj{marker,4}(m0:end1,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
                
                rt = plot3((store_Traj{marker,2}(m0:end1,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
                hold on
                
            end   
            
           
            
       end
    
        rt = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 1 0 0.5]);
    xl = xlabel('$x \,[$m$]$','Interpreter','latex');
%      xl.Position = [0.0112 -2.3821 -0.0445+0.002];
        xl.Position=[-0.1256 -5.1368 -0.0523];
    xl.HorizontalAlignment = 'center';

    yl =  ylabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');
%       yl.Position = [-0.5850 0.3355 -0.0449+0.001];
yl.Position = [-0.7871 -2.4698 -0.0546];
%     yl.Position = [-0.6460 1 -0.2637];

    yl.HorizontalAlignment = 'center';

    zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    zl.HorizontalAlignment = 'center';
    zl.Position = [-0.5567 1.5177 -0.0318];
    zticks([-0.05 -0.02 0])
%        rt = plot3(0,0,-0.1,'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
 title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$, $\alpha =',num2str(epsilon*round(t_ssm(1,1),4),'%4.4f'),'$');
 title(title_string,'Interpreter','latex')
       
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%        legend([h ft rt unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Stable solutions','Interpreter','latex','location','northwest')
    l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(\alpha))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.1474 0.7127 0.4769 0.2188];
       %       legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%      axis([-0.7,0.7,-1.5,1.5,-0.16,0.02])
%     view(-79,7)
set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')


    grid on
    box on
%            daspect([1,1,1])
   axis([-0.5,0.5,-1.5,1.5,-0.06,0])
 view(-10,10)
% axis([-0.6,0.6,-1.8,1.8,-0.1,0])
% view(-39,43)
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

     
        

     end 
   


end

for j = 120:20:1667
      clf
              plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-','LineWidth',2,'Color',[1 0 0 0.5]);
       hold on 
 
%          plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-.','LineWidth',1,'Color',[1 0 0 0.5]);
%        hold on 
         plotdofs = [2 4 3]; 
     m0 = j;
     x = linspace(-1,1,nssm);
     y = linspace(-0.5,0.5,nssm);
     [X,Y]=meshgrid(x,y);
    
     Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{1,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
             unst = plot3(store_Traj{1,2}(j:j+tl,plotdofs(1)),store_Traj{1,2}(j:j+tl,plotdofs(2)),store_Traj{1,2}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+tl,plotdofs(1)),store_Traj{1,2}(j+tl,plotdofs(2)),store_Traj{1,2}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       
        rt = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 1 0 0.5]);
        ft = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 0 0 0.5]);
    xl = xlabel('$x \,[$m$]$','Interpreter','latex');
%      xl.Position = [0.0112 -2.3821 -0.0445+0.002];
        xl.Position=[-0.1256 -5.1368 -0.0523];
    xl.HorizontalAlignment = 'center';

    yl =  ylabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');
%       yl.Position = [-0.5850 0.3355 -0.0449+0.001];
yl.Position = [-0.7871 -2.4698 -0.0546];
%     yl.Position = [-0.6460 1 -0.2637];

    yl.HorizontalAlignment = 'center';

    zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    zl.HorizontalAlignment = 'center';
    zl.Position = [-0.5567 1.5177 -0.0318];
    zticks([-0.05 -0.02 0])
%        rt = plot3(0,0,-0.1,'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
 title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$, $\alpha =',num2str(epsilon*round(t_ssm(1,1),4),'%4.4f'),'$');
 title(title_string,'Interpreter','latex')
       
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%        legend([h ft rt unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Stable solutions','Interpreter','latex','location','northwest')
    l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(\alpha))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.1474 0.7127 0.4769 0.2188];
       %       legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%      axis([-0.7,0.7,-1.5,1.5,-0.16,0.02])
%     view(-79,7)
set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')


    grid on
    box on
%            daspect([1,1,1])
   axis([-0.5,0.5,-1.5,1.5,-0.06,0])
 view(-10,10)
% axis([-0.6,0.6,-1.8,1.8,-0.1,0])
% view(-39,43)
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;
         
end    

for marker = 10:12 %to(4)+10:to(4)+10%1:spand:17000 % 1:100:70000%:1200%(max(size(t_sol))-tl)

     
% subplot(2, 2, [1 3])
   
     

       
%       plot3(store_Traj{1,3}(1:50:end,plotdofs(1)),store_Traj{1,3}(1:50:end,plotdofs(2)),store_Traj{1,3}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(1:50:end,plotdofs(1)),store_Traj{1,4}(1:50:end,plotdofs(2)),store_Traj{1,4}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
     for j = 1:40
          clf
         plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-','LineWidth',2,'Color',[1 0 0 0.5]);
        hold on 
         plotdofs = [2 4 3]; 
     m0 = j;
     x = linspace(-1,1,nssm);
     y = linspace(-0.5,0.5,nssm);
     [X,Y]=meshgrid(x,y);
     
    if m0+tl <=40
                  Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
             unst = plot3(store_Traj{1,2}(j+to2(marker-8-1),plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to2(marker-8-1),plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  

%               unst = plot3(store_Traj{1,2}(j+to2(marker-8-1)-1:j+to2(marker-8-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1)-1:j+to2(marker-8-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1)-1:j+to2(marker-8-1)-1+tl,plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
% %       hold on 
% %       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       
% %       hold on
%        plot3(store_Traj{1,2}(j+to2(marker-8-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1)-1+tl,plotdofs(3)),'.','MarkerSize',30,'Color','red');
%       hold on 
% %       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
% %       
      elseif m0<=40
                 Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
          unst = plot3(store_Traj{1,2}(j+to2(marker-8-1),plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to2(marker-8-1),plotdofs(1)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(2)),store_Traj{1,2}(j+to2(marker-8-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  
    end

       for nTraj = 1:16
             
             if m0+tl <=40
                 
                 ft = plot3(store_Traj{marker,4}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
            rt = plot3((store_Traj{marker,2}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
             plot3((store_Traj{marker,2}(m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(3)+(nTraj-1)*4)),'.','MarkerSize',30,'Color','green');
            hold on
             plot3(store_Traj{marker,4}(m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','black');
            hold on
            elseif m0<=40
                
                end1=40;
                ft = plot3(store_Traj{marker,4}(m0:end1,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
                
                rt = plot3((store_Traj{marker,2}(m0:end1,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
                hold on
                
            end   
            
           
            
       end
%       
        rt = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 1 0 0.5]);
    xl = xlabel('$x \,[$m$]$','Interpreter','latex');
%      xl.Position = [0.0112 -2.3821 -0.0445+0.002];
        xl.Position=[-0.1256 -5.1368 -0.0523];
    xl.HorizontalAlignment = 'center';

    yl =  ylabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');
%       yl.Position = [-0.5850 0.3355 -0.0449+0.001];
yl.Position = [-0.7871 -2.4698 -0.0546];
%     yl.Position = [-0.6460 1 -0.2637];

    yl.HorizontalAlignment = 'center';

    zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    zl.HorizontalAlignment = 'center';
    zl.Position = [-0.5567 1.5177 -0.0318];
    zticks([-0.05 -0.02 0])
%        rt = plot3(0,0,-0.1,'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
 title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$, $\alpha =',num2str(epsilon*round(t_ssm(1,1),4),'%4.4f'),'$');
 title(title_string,'Interpreter','latex')
       
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%        legend([h ft rt unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Stable solutions','Interpreter','latex','location','northwest')
    l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(\alpha))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.1474 0.7127 0.4769 0.2188];
       %       legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%      axis([-0.7,0.7,-1.5,1.5,-0.16,0.02])
%     view(-79,7)
set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')


    grid on
    box on
%            daspect([1,1,1])
   axis([-0.5,0.5,-1.5,1.5,-0.06,0])
 view(-10,10)
% axis([-0.6,0.6,-1.8,1.8,-0.1,0])
% view(-39,43)
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

     
        

     end 
   


end

for j = 1667+120:20:2500
      clf
              plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-','LineWidth',2,'Color',[1 0 0 0.5]);
       hold on 
 
%          plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-.','LineWidth',1,'Color',[1 0 0 0.5]);
%        hold on 
         plotdofs = [2 4 3]; 
     m0 = j;
     x = linspace(-1,1,nssm);
     y = linspace(-0.5,0.5,nssm);
     [X,Y]=meshgrid(x,y);
    
     Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{1,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
             unst = plot3(store_Traj{1,2}(j:j+tl,plotdofs(1)),store_Traj{1,2}(j:j+tl,plotdofs(2)),store_Traj{1,2}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+tl,plotdofs(1)),store_Traj{1,2}(j+tl,plotdofs(2)),store_Traj{1,2}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       
        rt = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 1 0 0.5]);
        ft = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 0 0 0.5]);
    xl = xlabel('$x \,[$m$]$','Interpreter','latex');
%      xl.Position = [0.0112 -2.3821 -0.0445+0.002];
        xl.Position=[-0.1256 -5.1368 -0.0523];
    xl.HorizontalAlignment = 'center';

    yl =  ylabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');
%       yl.Position = [-0.5850 0.3355 -0.0449+0.001];
yl.Position = [-0.7871 -2.4698 -0.0546];
%     yl.Position = [-0.6460 1 -0.2637];

    yl.HorizontalAlignment = 'center';

    zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    zl.HorizontalAlignment = 'center';
    zl.Position = [-0.5567 1.5177 -0.0318];
    zticks([-0.05 -0.02 0])
%        rt = plot3(0,0,-0.1,'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
 title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$, $\alpha =',num2str(epsilon*round(t_ssm(1,1),4),'%4.4f'),'$');
 title(title_string,'Interpreter','latex')
       
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%        legend([h ft rt unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Stable solutions','Interpreter','latex','location','northwest')
    l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(\alpha))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.1474 0.7127 0.4769 0.2188];
       %       legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%      axis([-0.7,0.7,-1.5,1.5,-0.16,0.02])
%     view(-79,7)
set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')


    grid on
    box on
%            daspect([1,1,1])
   axis([-0.5,0.5,-1.5,1.5,-0.06,0])
 view(-10,10)
% axis([-0.6,0.6,-1.8,1.8,-0.1,0])
% view(-39,43)
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;
         
end    



for marker = 18:20 %to(4)+10:to(4)+10%1:spand:17000 % 1:100:70000%:1200%(max(size(t_sol))-tl)

     
% subplot(2, 2, [1 3])
   
     

       
%       plot3(store_Traj{1,3}(1:50:end,plotdofs(1)),store_Traj{1,3}(1:50:end,plotdofs(2)),store_Traj{1,3}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(1:50:end,plotdofs(1)),store_Traj{1,4}(1:50:end,plotdofs(2)),store_Traj{1,4}(1:50:end,plotdofs(3)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
     for j = 1:40
         clf
         plot3(store_Traj{1,2}(1:30:end,plotdofs(1)),store_Traj{1,2}(1:30:end,plotdofs(2)),store_Traj{1,2}(1:30:end,plotdofs(3)),'-','LineWidth',2,'Color',[1 0 0 0.5]);
       hold on 
         plotdofs = [2 4 3]; 
     m0 = j;
     x = linspace(-1,1,nssm);
     y = linspace(-0.5,0.5,nssm);
     [X,Y]=meshgrid(x,y);
     
    if m0+tl <=40
                  Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
                    unst = plot3(store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  

%                    unst = plot3(store_Traj{1,2}(j+to3(marker-2*8-1)-1:j+to3(marker-2*8-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1)-1:j+to3(marker-2*8-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1)-1:j+to3(marker-2*8-1)-1+tl,plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
% %       hold on 
% %       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
% %       
% %       hold on
%        plot3(store_Traj{1,2}(j+to3(marker-2*8-1)-1+tl,plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1)-1+tl,plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1)-1+tl,plotdofs(3)),'.','MarkerSize',30,'Color','red');
%       hold on 
% %       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
% %       hold on 
% %       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
% %       

      elseif m0<=40
                 Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*store_Traj{marker,1}(j);

     [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);

  
     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);


     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 
               unst = plot3(store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(3)),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%       hold on 
%       st = plot3(store_Traj{1,3}(j:j+tl,plotdofs(1)),store_Traj{1,3}(j:j+tl,plotdofs(2)),store_Traj{1,3}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3(store_Traj{1,4}(j:j+tl,plotdofs(1)),store_Traj{1,4}(j:j+tl,plotdofs(2)),store_Traj{1,4}(j:j+tl,plotdofs(3)),'-','LineWidth',3,'Color',[0 0 1 0.5]);
%       
%       hold on
       plot3(store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(1)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(2)),store_Traj{1,2}(j+to3(marker-2*8-1),plotdofs(3)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3(store_Traj{1,3}(j+tl,plotdofs(1)),store_Traj{1,3}(j+tl,plotdofs(2)),store_Traj{1,3}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3(store_Traj{1,4}(j+tl,plotdofs(1)),store_Traj{1,4}(j+tl,plotdofs(2)),store_Traj{1,4}(j+tl,plotdofs(3)),'.','MarkerSize',30,'Color',[0,0,1]);
%  

    end

       for nTraj = 1:16
             
             if m0+tl <=40
                 
                 ft = plot3(store_Traj{marker,4}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
            rt = plot3((store_Traj{marker,2}(m0:m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:m0+tl,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
             plot3((store_Traj{marker,2}(m0+tl,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0+tl,plotdofs(3)+(nTraj-1)*4)),'.','MarkerSize',30,'Color','green');
            hold on
             plot3(store_Traj{marker,4}(m0+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','black');
            hold on
            elseif m0<=40
                
                end1=40;
                ft = plot3(store_Traj{marker,4}(m0:end1,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,4}(m0:end1,plotdofs(3)+(nTraj-1)*4),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
                
                rt = plot3((store_Traj{marker,2}(m0:end1,plotdofs(1)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(2)+(nTraj-1)*4)),(store_Traj{marker,2}(m0:end1,plotdofs(3)+(nTraj-1)*4)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
                hold on
                
            end   
            
           
            
       end
        rt = plot3(0,0,-0.01,'-','LineWidth',3,'Color',[0 1 0 0.5]);
    xl = xlabel('$x \,[$m$]$','Interpreter','latex');
%      xl.Position = [0.0112 -2.3821 -0.0445+0.002];
        xl.Position=[-0.1256 -5.1368 -0.0523];
    xl.HorizontalAlignment = 'center';

    yl =  ylabel('$\dot{x} \,[$m/s$]$','Interpreter','latex');
%       yl.Position = [-0.5850 0.3355 -0.0449+0.001];
yl.Position = [-0.7871 -2.4698 -0.0546];
%     yl.Position = [-0.6460 1 -0.2637];

    yl.HorizontalAlignment = 'center';

    zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
    zl.HorizontalAlignment = 'center';
    zl.Position = [-0.5567 1.5177 -0.0318];
    zticks([-0.05 -0.02 0])
%        rt = plot3(0,0,-0.1,'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
 title_string = strcat('$\epsilon = ',num2str(epsilon),'$, expansion order $N = 3$, $\alpha =',num2str(epsilon*round(t_ssm(1,1),4),'%4.4f'),'$');
 title(title_string,'Interpreter','latex')
       
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%        legend([h ft rt unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Stable solutions','Interpreter','latex','location','northwest')
    l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(\alpha))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.1474 0.7127 0.4769 0.2188];
       %       legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%      axis([-0.7,0.7,-1.5,1.5,-0.16,0.02])
%     view(-79,7)
set(gcf,'color','white')
set(hFig, 'Units' , 'Inches' );
pos = get(hFig, 'Position' );
set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
set(gcf,'Renderer','painters')


    grid on
    box on
%            daspect([1,1,1])
   axis([-0.5,0.5,-1.5,1.5,-0.06,0])
 view(-10,10)
% axis([-0.6,0.6,-1.8,1.8,-0.1,0])
% view(-39,43)
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

     
        

     end 
   


end

myWriter = VideoWriter('Bumpy_rail_SSM_adiabatic','MPEG-4');
myWriter.FrameRate = 10;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);

% set(hFig, 'Units' , 'Inches' );
% pos = get(hFig, 'Position' );
% set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(hFig, 'ssm_rail_slow_snap1.pdf' , '-dpdf')

% set(fig, 'Units' , 'Inches' );
% pos = get(fig, 'Position' );
% set(fig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(fig, 'forcing_slow.pdf' , '-dpdf')
% 

% plot3(SP_Traj1(plotdofs(1),j:j+tl),SP_Traj1(plotdofs(2),j:j+tl),SP_Traj1(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color',[0 0 0 0.5]);
%      hold on
%      plot3(shift_traj(plotdofs(1),j:j+tl),shift_traj(plotdofs(2),j:j+tl),shift_traj(plotdofs(3),j:j+tl),'-','LineWidth',3,'Color',[1 0 0 0.5]);
%      hold on
%      plot3(SolutionA(plotdofs(1),j:j+tl),SolutionA(plotdofs(2),j:j+tl),SolutionA(plotdofs(3),j:j+tl),'-','LineWidth',3,'color',[0,0,1,0.8])
%      hold on
%      plot3(shift_traj(plotdofs(1),j+tl),shift_traj(plotdofs(2),j+tl),shift_traj(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','red');
%      hold on
%      plot3(SP_Traj1(plotdofs(1),j+tl),SP_Traj1(plotdofs(2),j+tl),SP_Traj1(plotdofs(3),j+tl),'.','MarkerSize',30,'Color','black');
%      hold on
%      plot3(SolutionA(plotdofs(1),j+tl),SolutionA(plotdofs(2),j+tl),SolutionA(plotdofs(3),j+tl),'.','MarkerSize',30,'Color',[0,0,1]);
%      hold on