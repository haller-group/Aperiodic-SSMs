function IC = getStaticResponse(K, M, F, loadvector, plotMesh, PlotFieldonDefMesh)

nTraj = size(loadvector, 2);
n = size(K, 1);
w0 = -K\loadvector(:,1); % linear initial guess
IC = zeros(2*n,nTraj);
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 100000*n,'MaxIterations', 100000*n, 'Display', 'off');

for iLoad = 1:nTraj
    w0 = -K\loadvector(:,nTraj);
    f_eq = @(w)([zeros(n) M]*F(0,[w; zeros(n,1)]) + loadvector(:,iLoad));
    tic
    [w0, ~, exitflag, output] = fsolve(f_eq, w0, options);
    toc
    if exitflag <= 0
        disp(output);
        error('Warning: No solution found for loading configuration')
    end
    IC(:,iLoad) = [w0; zeros(n,1)];
    
%     PlotFieldonDefMesh(w0,0.2)
%     hold on 
end
% if plotMesh
%     figure; 
%     PlotFieldonDefMesh(w0,200)
% end