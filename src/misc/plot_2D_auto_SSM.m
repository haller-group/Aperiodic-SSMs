function [zdof1,zdof2,zdof3]=plot_2D_auto_SSM(coord_change_vec2,coord_change_vec,W,rhosamp,plotdofs,varargin)
% PLOT_2D_AUTO_SSM This function returns the plot of 2 dimensional
% autonomous SSM. We first generate the grids of parameterization
% coordinates based on rhosamp X thetasamp(0:2*pi/128:2*pi). Then we map
% these grids to the full system based on the expansion W of SSM at
% plotdofs. Note that plotdofs should have 3 dofs.

assert(numel(plotdofs)==3,'the number of dofs is not three');
% generate grids
[RHO,THETA] = meshgrid(rhosamp,0:2*pi/128:2*pi);
zdof1 = zeros(size(RHO));
zdof2 = zeros(size(RHO));
zdof3 = zeros(size(RHO));
for k=1:129
    pk = RHO(k,:).*exp(1i*THETA(k,:));
    zk = reduced_to_full_traj([],[pk;conj(pk)],W);
    zk = zk(plotdofs,:);
%     zdof1(k,:) = zk(1,:);
%     zdof2(k,:) = zk(2,:);
%     zdof3(k,:) = zk(3,:);
       zdof1(k,:) = zk(1,:)+coord_change_vec(plotdofs(1));
    zdof2(k,:) = zk(2,:)+coord_change_vec(plotdofs(2));
    zdof3(k,:) = zk(3,:)+coord_change_vec(plotdofs(3));
%     zdof1(k,:) = zk(2,:)+coord_change_vec(plotdofs(2));
%     zdof2(k,:) = zk(2,:)+coord_change_vec(plotdofs(1));
%     zdof3(k,:) = zk(3,:)+coord_change_vec(plotdofs(3));

end


h = surf(zdof1,zdof2,zdof3);

h.EdgeColor = 'k';
h.FaceColor = [0.7 0.7 0.7];
h.FaceAlpha = 0.3;

view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',30);
if isempty(varargin)
    xlabel(['$z_{\mathrm{',num2str(plotdofs(1)),'}}$'],'interpreter','latex');
    ylabel(['$z_{\mathrm{',num2str(plotdofs(2)),'}}$'],'interpreter','latex');
    zlabel(['$z_{\mathrm{',num2str(plotdofs(3)),'}}$'],'interpreter','latex');
else
    xlabel(varargin{1}{1},'interpreter','latex');
    ylabel(varargin{1}{2},'interpreter','latex');
    zlabel(varargin{1}{3},'interpreter','latex');
%     title('SSMs in the phase space anchored to a target','interpreter','latex');
end    

end


