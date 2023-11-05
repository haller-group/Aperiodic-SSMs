function RF = grid_SSM_interp1_nc(KT,X,Ngrid)

RA = cell(max(size(KT{1,1})),1);
for k = 1:max(size(KT{1,1}))
RA{k,1} = [];
end
for i = 1:Ngrid
    for j = 1:max(size(KT{1,i}))
        RA{j,1} =[RA{j,1}, KT{1,i}(:) ]; 
    end    
end    

% RF=0;


% fit interpolant spline
RF = cell(max(size(KT{1,1})),1);
for ind = 1:max(size(RA))
   if  ~isempty(RA{ind,1})
        PF = RA{ind,1};
        FrF = cell(max(size(PF(:,1))),1);
        FiF = cell(max(size(PF(:,1))),1);
        for indj = 1:max(size(PF(:,1)))
            L = PF(indj,:);
            
            Lr = full(L);
            Fr = griddedInterpolant(X,real(Lr),'spline','spline');
            Fi = griddedInterpolant(X,imag(Lr),'spline','spline');
            FrF(indj,1) = {Fr};
            FiF(indj,1) = {Fi};
        end    
       
        RF{ind,1} = {FrF};
   else
        RF{ind,1} = [];
   end    
end

end