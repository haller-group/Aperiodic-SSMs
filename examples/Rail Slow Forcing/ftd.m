function [dXdt, Xtrunc,t_trunc] = ftd(X, t) % finite time difference
    ind = 5:size(X,2)-4; Xtrunc = X(:,ind); t_trunc = t(ind);
    dX = 4/5*(X(:,ind+1)-X(:,ind-1)) - 1/5*(X(:,ind+2)-X(:,ind-2)) + ...
        4/105*(X(:,ind+3)-X(:,ind-3)) - 1/280*(X(:,ind+4)-X(:,ind-4));
    dXdt = dX./(t(ind+1)-t(ind));
end
