hFig = figure('DefaultAxesFontSize',18);                       % Bring up new figure
% imshow('board.tif','Border','tight') % The axes will fill up the entire figure as much as possible without changing aspect ratio.
% hFig.WindowState = 'maximized';
% pause(0.5)
delta = 0.5;
tl =  10;
% tl = 1000
clear movieVector

nssm = 20;
m1 = 1+3;
m2 = 1;
m3 = 1;
spand = 1;
ind = 1;
for j = to(3):to(3) %1:20 %1:20 %1:spand:30 % 1:100:70000%:1200%(max(size(t_sol))-tl)
i = j+0;
     clf
% subplot(2, 2, [1 3])
    
     plotdofs = [1 3 4]; 
   
     x = linspace(-0.6,0.6,nssm);
     y = linspace(-1.4,1.4,nssm);
     [X,Y]=meshgrid(x,y);
    
     Z = X+1i*Y;
     Z = [X(:),Y(:)];
     t_ssm = ones(1,nssm*nssm)*store_Traj{1,1}(j+tl);

%      [ZZ,modalZ,Net_Solz] = compute_SSM_phy(XI_0U,XI_1U,XI_2U,XI_3U,t_ssm,Z(:).',epsilon,SSM_Coeff_A_2,SSM_Coeff_A_1,Valpha,V,A,Force_Lorenz,Dalpha);
%      [ZZ,modalZ] = compute_SSM_phy(XI_hyperbolic_solution,t_ssm,V_new,VI_new,Z.',h030,h300,h120,h210,h111,h201,h021,h012,h102,h500,h050,h410,h140,h230,h320,epsilon);

  [ZZ,modalZ] = compute_SSM_phy(XI_hyperbolic_solution,t_ssm,V_new,VI_new,Z.',h030,h300,h120,h210,h111,h201,h021,h012,h102,epsilon);

%      Z1 = reshape(real(modalZ(plotdofs(1),:)),nssm,nssm);
%      Z2 = reshape(real(modalZ(plotdofs(2),:)),nssm,nssm);
%      Z3 = reshape(imag(modalZ(plotdofs(3),:)),nssm,nssm);

     Z1 = reshape(ZZ(plotdofs(1),:),nssm,nssm);
     Z2 = reshape(ZZ(plotdofs(2),:),nssm,nssm);
     Z3 = reshape(ZZ(plotdofs(3),:),nssm,nssm);



     h = surf(Z1,Z2,Z3)
     h.EdgeColor = [0 0 0];
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     h.EdgeAlpha = 0.3;
     hold on 

%         plot3((store_Traj{1,2}(plotdofs(1),1:2000)),(store_Traj{1,2}(plotdofs(2),1:2000)),(store_Traj{1,2}(plotdofs(3),1:2000)),'-.','LineWidth',1,'Color',[1 0 0 0.5]);
%        hold on 
%       plot3((store_Traj{1,3}(plotdofs(1),1500:2000)),(store_Traj{1,3}(plotdofs(2),1500:2000)),(store_Traj{1,3}(plotdofs(3),1500:2000)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%       hold on 
%       plot3((store_Traj{1,4}(plotdofs(1),1500:2000)),(store_Traj{1,4}(plotdofs(2),1500:2000)),(store_Traj{1,4}(plotdofs(3),1500:2000)),'-.','LineWidth',1,'Color',[0 0 1 0.5]);
%        hold on 
%      if (j+tl >= to(1))
%          marker = 2;
%        for nTraj = 1:16
%           
%              ft = plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,j:j+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,j:j+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,j:j+tl)),'-','LineWidth',3,'Color',[0 0 0 0.5]);
%             hold on 
%           
%             if j+tl <=200
%             rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,j:j+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,j:j+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,j:j+tl)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%             hold on 
%              plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,j+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,j+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,j+tl)),'.','MarkerSize',30,'Color','green');
%             hold on
%             elseif j<=200
%                 end1=200;
%                 rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,j:end1)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,j:end1)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,j:end1)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%                 hold on
%             end    
%             
%            
%             plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,j+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,j+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,j+tl)),'.','MarkerSize',30,'Color','black');
%             hold on
% %             plot3(store_Traj{marker,2}(j+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','green');
% %             hold on
%        end  
%      end  
%      if   (j+tl >= to(2))
%            marker = 3;
% 
%             for nTraj = 1:16
%           
%              ft = plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m1:m1+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m1:m1+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m1:m1+tl)),'-','LineWidth',3,'Color',[0 0 0 0.5]);
% 
%             if m1+tl <=26
%             rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m1:m1+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m1:m1+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m1:m1+tl)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%             hold on 
%              plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m1+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m1+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m1+tl)),'.','MarkerSize',30,'Color','green');
%             hold on
%             elseif m1<=26
%                 end1=26;
%                 rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m1:end1)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m1:end1)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m1:end1)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%                 hold on
%             end 
% 
%             plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m1+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m1+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m1+tl)),'.','MarkerSize',30,'Color','black');
%             hold on
% 
% %             plot3(store_Traj{marker,2}(j+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','green');
% %             hold on
%            end 
%       
%        m1 = m1 + spand;
%      end   
     if     (j+tl >= to(3))
          marker = 4;
          for nTraj = 1:16
          
             ft = plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m2:m2+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m2:m2+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m2:m2+tl)),'-','LineWidth',3,'Color',[0 0 0 0.5]);
            hold on 
             if m2+tl <=200
            rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m2:m2+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m2:m2+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m2:m2+tl)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
            hold on 
             plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m2+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m2+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m2+tl)),'.','MarkerSize',30,'Color','green');
            hold on
            elseif m2<=200
                end1=200;
                rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m2:end1)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m2:end1)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m2:end1)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
                hold on
            end

            plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m2+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m2+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m2+tl)),'.','MarkerSize',30,'Color','black');
            hold on

%             plot3(store_Traj{marker,2}(j+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','green');
%             hold on
           end 
        
       m2 = m2 + spand;
     end
%      if (j+tl >= to(4))
%           marker = 5;
%        for nTraj = 1:16
%           
%             ft = plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m3:m3+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m3:m3+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m3:m3+tl)),'-','LineWidth',3,'Color',[0 0 0 0.5]);
%             hold on 
%              if m3+tl <=200
%             rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m3:m3+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m3:m3+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m3:m3+tl)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%             hold on 
%              plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m3+tl)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m3+tl)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m3+tl)),'.','MarkerSize',30,'Color','green');
%             hold on
%             elseif m3<=200
%                 end1=200;
%                 rt = plot3((store_Traj{marker,2}(plotdofs(1)+(nTraj-1)*4,m3:end1)),(store_Traj{marker,2}(plotdofs(2)+(nTraj-1)*4,m3:end1)),(store_Traj{marker,2}(plotdofs(3)+(nTraj-1)*4,m3:end1)),'-','LineWidth',3,'Color',[0 1 0 0.5]);
%                 hold on
%             end
% 
%             plot3((store_Traj{marker,4}(plotdofs(1)+(nTraj-1)*4,m3+tl)),(store_Traj{marker,4}(plotdofs(2)+(nTraj-1)*4,m3+tl)),(store_Traj{marker,4}(plotdofs(3)+(nTraj-1)*4,m3+tl)),'.','MarkerSize',30,'Color','black');
%             hold on
% %             plot3(store_Traj{marker,2}(j+tl,plotdofs(1)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(2)+(nTraj-1)*4),store_Traj{marker,2}(j+tl,plotdofs(3)+(nTraj-1)*4),'.','MarkerSize',30,'Color','green');
% %             hold on
%         end  
%         
%        m3 = m3 + spand;
% 
%      end 
     hold on
      unst = plot3((store_Traj{1,2}(plotdofs(1),j:j+tl)),(store_Traj{1,2}(plotdofs(2),j:j+tl)),(store_Traj{1,2}(plotdofs(3),j:j+tl)),'-','LineWidth',3,'Color',[1 0 0]);
      hold on 
%       st = plot3((store_Traj{1,3}(plotdofs(1),j:j+tl)),(store_Traj{1,3}(plotdofs(2),j:j+tl)),(store_Traj{1,3}(plotdofs(3),j:j+tl)),'-','LineWidth',3,'Color',[0 0 1]);
%       hold on 
%       plot3((store_Traj{1,4}(plotdofs(1),j:j+tl)),(store_Traj{1,4}(plotdofs(2),j:j+tl)),(store_Traj{1,4}(plotdofs(3),j:j+tl)),'-','LineWidth',3,'Color',[0 0 1]);
%       
%       hold on
       plot3((store_Traj{1,2}(plotdofs(1),j+tl)),(store_Traj{1,2}(plotdofs(2),j+tl)),(store_Traj{1,2}(plotdofs(3),j+tl)),'.','MarkerSize',30,'Color','red');
      hold on 
%       plot3((store_Traj{1,3}(plotdofs(1),j+tl)),(store_Traj{1,3}(plotdofs(2),j+tl)),(store_Traj{1,3}(plotdofs(3),j+tl)),'.','MarkerSize',30,'Color',[0,0,1]);
%       hold on 
%       plot3((store_Traj{1,4}(plotdofs(1),j+tl)),(store_Traj{1,4}(plotdofs(2),j+tl)),(store_Traj{1,4}(plotdofs(3),j+tl)),'.','MarkerSize',30,'Color',[0,0,1]);
      
      
    xl =xlabel('$x\,[$m$]$','Interpreter','latex');

%     xl.Position = [0.0793 -2.4921+0.2 -0.1372]
%      xl.Position = [0.0793 -2.2921 -0.0072];

      xl.Position = [0.0485 -2.3265 1.0424]

% xl.Position = [0.1111 -1.0188+0.05 -1.9472+0.05]
    xl.HorizontalAlignment = 'center';
     yl =ylabel('$\dot{x}\,[$m/s$]$','Interpreter','latex');
%      yl.Position = [-0.5002 0.5250 -0.1382];
%       yl.Position = [-0.5002 0.5250 -0.0082];
       yl.Position = [-0.4539 +0.2 1.0421];


    yl.HorizontalAlignment = 'center';
     zl = zlabel('$\dot{x}_c \,[$m/s$]$','Interpreter','latex');
     title_string = strcat('max$(|F_{ext}|) = ',num2str(epsilon*(M+m)),'$ [$\mathbf{N}$], expansion order $N = 3$, $t =',num2str(round(t_ssm(1,1),2),'%4.2f'),'$');
       title(title_string,'Interpreter','latex')
%          legend('SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','$O(\epsilon^3)$ anchor trajectory','Full order model','Reduced order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','Interpreter','latex','location','northwest')
%   rt = plot3(0,0,0,'-','LineWidth',3,'Color',[0 1 0 0.5]);
%             hold on 
% l = legend([h ft rt unst st],'SSM - $W_{\epsilon}(x^{*}(t))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','Reduced order model','Anchor trajectory','Stable solutions','Interpreter','latex','location','northwest');
% l.Position = [0.1353 0.7847 0.2358 0.1311];
   l = legend([h ft rt unst],'Mixed mode SSM $W_{\epsilon}(x^{*}(t))$','Full order model','Reduced order model','Chaotic anchor trajectory','Interpreter','latex','location','northwest');
l.Position = [0.0411 0.7363 0.4284 0.1982];

%        legend([h ft unst st],'SSM - $W_{\epsilon}(x_{\epsilon}(\alpha))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 3)$','Full order model','$O(\epsilon^3)$ anchor trajectory [$x_{\epsilon}(\alpha = \epsilon t)$]','$O(\epsilon^3)$ stable solutions','Interpreter','latex','location','northwest')

%     string = strcat('\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),', \epsilon = 0.008');
%     if (round(epsilon*t_sol(j+tl),2)) >=1
%      string = strcat('{\color{blue}\alpha = ',num2str(round(epsilon*t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.008');
%      title(string);
%     else
%      title(string);
%     end   
%        axis([-0.2,0.2,-1,1,-2*10^-3,18*10^-3])
%    view(-34,6)


    grid on
    box on
%             daspect([1,1,1])
%    axis([-0.4,0.4,-2,2,-0.1,0.24])
set(gcf,'color','white')
    figssm = gcf;

movieVector(ind) = getframe(hFig);
    
    ind = ind +1;

end

% myWriter = VideoWriter('SSM_Rail_solutions_weak_slow','MPEG-4');
% myWriter.FrameRate = 20;
% 
% open(myWriter);
% writeVideo(myWriter,movieVector);
% close(myWriter);

% 
% 
% for indj = 1:732
% writeVideo(myWriter,movieVector(indj));
% end
% 
% close(myWriter);

% set(hFig, 'Units' , 'Inches' );
% pos = get(hFig, 'Position' );
% set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(hFig, 'ssm_rail_weak_snap2.pdf' , '-dpdf')


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