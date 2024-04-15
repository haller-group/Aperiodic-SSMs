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

for j =1:1000% 150:2700 %2550:2550 %1000:1000 %150:150 %150:2700%2550 %1:20:5000%(max(size(t_sol))-tl)
i = j+0;

     clf
     rhosamp = linspace(0,0.5,201);
     plotdofs = [1 2 4]; 
     [RHO,Theta]=meshgrid(linspace(0,1,201),linspace(0,2*pi,201));
     x = linspace(-2,2,nssm);
     [X,Y]=meshgrid(x);
     Z = RHO.*exp(1i*Theta);
     Z = X+1i*Y;
     t_ssm = ones(1,nssm*nssm)*t_sol(j+tl);
     [ZZ,modalZ] = compute_SSM_phy(xi_full_eval5,xi_full_eval3,xi_full_eval,t_ssm,l1,l1,gamma,V,inv(V),Z(:).',h030,h300,h210,h120,H_Coeff1,f030,f300,f210,f120,F_Coeff1,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,H_Coeff2,F_Coeff2,epsilon);
    
     Z1 = reshape(ZZ(1,:),nssm,nssm);
     Z2 = reshape(ZZ(2,:),nssm,nssm);
     Z3 = reshape(ZZ(3,:),nssm,nssm);


     starto = 1;
     Z1 = reshape(real(modalZ(1,:)),nssm,nssm);
     Z2 = reshape(imag(modalZ(1,:)),nssm,nssm);
     Z3 = reshape(real(modalZ(3,:)),nssm,nssm);

     h = surf(Z1,Z2,Z3);
     h.EdgeColor = 'k';
     h.FaceColor = [.7 .7 .7];
     h.FaceAlpha = 0.3;
     hold on 

     starto = 1;
     
     hold on
     plot3(SP_TrajZ(plotdofs(1),j:i+tl),SP_TrajZ(plotdofs(2),j:i+tl),SP_TrajZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'Color',[0 0 0 0.5]);
     hold on
     plot3(shift_trajZ(plotdofs(1),j:i+tl),shift_trajZ(plotdofs(2),j:i+tl),shift_trajZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'Color',[1,0,0,0.5]);
     hold on
     plot3(SolutionZ(plotdofs(1),j:i+tl),SolutionZ(plotdofs(2),j:i+tl),SolutionZ(plotdofs(3),j:i+tl),'-','LineWidth',3,'color',[0 0 1 0.5])
     hold on
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

     
      legend('SSM - $W_{\epsilon}(x_{\epsilon}(t))$ $O(\epsilon^{k_1} \mathbf{u}^{\mathbf{k}_2}, k_1 + |k_2| = 5)$','Full order model','Reduced order model','Chaotic anchor trajectory ($x_{\epsilon}(t)$)','Interpreter','latex','location','northeast')

    string = strcat('$t = $',num2str(round(t_sol(j+tl),2),'%4.2f' ),', max$(|F_{ext}|) = 3$ [$\mathbf{N}$]');
    if (j+tl) >=3000
     string = strcat('{\color{blue}t = ',num2str(round(t_sol(j+tl),2),'%4.2f' ),',}\epsilon = 0.5');
     title(string);
    else
     title(string,'Interpreter','latex');
    end    
 view(-58,22)
    grid on
    box on 
     daspect([1,1,1])
    axis([-2.5,2.5,-2.5,2.5,-1.5,1.5])
%    axis([-1,1,-1,1,-1,1])
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


myWriter = VideoWriter('Shake_Cart_SSM_video_dis_fast','MPEG-4');
myWriter.FrameRate = 50;

open(myWriter);
% 
for indj = 1:1000
writeVideo(myWriter,movieVector(indj));
end

close(myWriter);

% set(hFig, 'Units' , 'Inches' );
% pos = get(hFig, 'Position' );
% set(hFig, 'PaperPositionMode' , 'Auto' , 'PaperUnits' , 'Inches' , 'PaperSize' ,[pos(3), pos(4)])
% set(gcf,'Renderer','painters')
% print(hFig, 'ssm_shaking_cart_snap3.pdf' , '-dpdf' , '-r300' )
