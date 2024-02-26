function [d] = GetDataInitial(fit,nLevel)
% generate new initial data values according to hetrodesic noise structure given a specified noise level
% returns generated initial data point 
    d = zeros(1,length(fit));
        for j = 1:length(fit)
            d(j) = fit(j) + normrnd(0,fit(j)*nLevel(1).^2);
            while d(j) < 0
            % resample if value is negative
                d(j) = fit(j) + normrnd(0,fit(j)*nLevel(1).^2);
            end
       
      
        end
end
