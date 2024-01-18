function [d] = GetDataInitial(fit,nLevel)
    d = zeros(1,length(fit));
        for j = 1:length(fit)
            d(j) = fit(j) + normrnd(0,fit(j)*nLevel(1).^2);
            while d(j) < 0
                d(j) = fit(j) + normrnd(0,fit(j)*nLevel(1).^2);
            end
       
      
        end
end