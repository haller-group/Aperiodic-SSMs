function FT = convert_fit_to_value_x(FF,xa,F0)
        FT = F0;
        for ind = 1:max(size(FF))
            if ~isempty(FF{ind,1})
                rf = FF{ind,1}{1,1};
                irf = FF{ind,1}{1,2};
                realf = [];
                imagf = [];
                for indj = 1:max(size(rf))
                    realf = [realf;rf{indj,1}(xa)];
                    imagf = [imagf;irf{indj,1}(xa)];
                end
                coeff_m = reshape(realf+1j*imagf,max(size(FT{1,ind}.coeffs(:,1))),max(size(FT{1,ind}.coeffs(1,:))));
                FT{1,ind}.coeffs = sparse(coeff_m);          
            end
        end    
end