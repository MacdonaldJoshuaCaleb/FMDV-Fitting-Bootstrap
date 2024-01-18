function vpfit = ViralGrowthforCI(ses,index)
%ses = readtable('ViremiaContact.csv');
%ses = table2array(ses(index,3:end))
tt = [0,2, 4, 6, 9,12,28];

function yy = paramfun1(p,t)
    x0 = p(1);
    r = p(2);
    f = @(t) x0.*exp(r.*t);
   yy = f(t);
end 

virus = ses;
    for k = 1:length(virus)
        if virus(k) >= .75
            I = k;
            break
        end
    end
    virus = virus(1:I);
    tts = tt(1:I);
    t = linspace(tts(1),tts(end));
    lb = [2.1e-7, 2.81];
    ub = [.01,4.15];
    params = readtable('FitParamsContact.csv');
   params = table2array(params(:,4:10));
   p0 = [params(index,2),params(index,1)];
   %p0 = (lb+ub)./2
    [vpfit resnorm] = lsqcurvefit(@paramfun1,p0,tts,virus,lb,ub);
    %[vpfit resnorm] = lsqcurvefit(@paramfun1,vpfit,tts,virus,lb,ub);
    %[vpfit resnorm] = lsqcurvefit(@paramfun1,vpfit,tts,virus,lb,ub);
%     fit = paramfun1(vpfit,t);
%     figure
%     plot(t,fit,'r',tts,virus,'r*')
%     Virusfits(j,:) = vpfit; 
%     save('Virusfits.mat','Virusfits')
end