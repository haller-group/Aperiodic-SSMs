function FT = convert_fit_to_value_x_na(FF,xa,F0)
        FT = F0;
        for ind = 1:max(size(FF))
            if ~isempty(FF{ind,1})
                rf = FF{ind,1}{1,1};
                
                realf = [];
                
                for indj = 1:max(size(rf))
                    realf = [realf;rf{indj,1}(xa)];
                    
                end
                coeff_m = reshape(realf,max(size(FT(:,1))),max(size(FT(1,:))));
                FT = coeff_m;  
        end    
end