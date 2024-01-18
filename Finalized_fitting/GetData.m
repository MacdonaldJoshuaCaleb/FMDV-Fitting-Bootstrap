function [d] = GetData(fit1,fit2,fit3,nLevel)
    d1 = zeros(1,length(fit1));
    d2 = zeros(1,length(fit1));
    d3 = zeros(1,length(fit1));
        for j = 1:length(fit1)
            d1(j) = fit1(j) + normrnd(0,fit1(j)*nLevel(1).^2);
            while d1(j) < 0
                d1(j) = fit1(j) + normrnd(0,fit1(j)*nLevel(1).^2);
            end
            %%%
            d2(j) = fit2(j) + normrnd(0,fit2(j)*nLevel(2).^2);
            while d2(j) < 0
                d2(j) = fit2(j) + normrnd(0,fit2(j)*nLevel(2).^2);
            end
            %%%
            d3(j) = fit3(j) + normrnd(0,fit3(j)*nLevel(3).^2);
            while d3(j) < 0
                d3(j) = fit3(j) + normrnd(0,fit3(j)*nLevel(3).^2);
            end
      
        end
        d = [d1, d2, d3];
end
