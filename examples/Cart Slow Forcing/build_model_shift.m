% generalised - not yet - only works till cubic nonlinearities
function [K_shift,fnl_shift] = build_model_shift(K,fnl,X_s)
n = length(K);

% new modifications
% Ideally create it with the length to modify

fnl_s = cell(1,length(fnl));
for iT = 1:length(fnl)
    tensor_dim = num2cell(repmat(n,1,iT+1));
    fnl_s{1,iT} = sptensor([tensor_dim{1,:}]);
end    

for iT = 1:length(fnl)
    Fi = fnl{iT};
    idx = find(Fi);
    idx_n = [];
    idx_l = [];
    P = [];
    if ~isempty(idx)
        for kT = 1:max(size(idx(:,1)))
        
            P = [P;idx(kT,2:end)];
            for i = 1:iT
                P = [P;circshift(idx(kT,2:end),i)];
            end
            P_full = [idx(kT,1)*ones(max(size(P(:,1))),1),P];
            idx_n = [idx_n;P_full];
            idx_l = [idx_l;repmat(idx(kT,:),max(size(P(:,1))),1)];
            P = [];

        end

    % Changes in All Parts
        idx_na=num2cell(idx_n);
        idx_la=num2cell(idx_l);
        for pl = 1:max(size(idx_na(:,1)))
            for loop_c = 1:iT
            
                n_set = [idx_na{pl,2:end-loop_c}];
                i_set = num2cell([idx_na{pl,1:1},idx_na{pl,end-loop_c+1:end}]);
            
                fnl_s{1,loop_c}(i_set{1,:}) = fnl_s{1,loop_c}(i_set{1,:}) +Fi(idx_la{pl,:})*prod(X_s(n_set));
            end    
        end 
    end
    

end

fnl_shift = cell(1,length(fnl));
for iT = 1:length(fnl)-1
fnl_shift{1,iT} = fnl{1,iT} + fnl_s{1,iT+1};
end 
fnl_shift{1,end} = fnl{1,end};


K_shift = K + ((spmatrix(fnl_s{1,1}))+(spmatrix(fnl_s{1,1})).')/2;

end

%         n_set = [idx_na{pl,2:end-1}];
%         K_s(idx_na{pl,1},idx_na{pl,end}) = K_s(idx_na{pl,1},idx_na{pl,end}) + Fi(idx_la{pl,:})*prod(X_s(n_set));
