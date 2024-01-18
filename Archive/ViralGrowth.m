function vpfit = ViralGrowth(ses)
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
        if virus(k) > 1
            I = k;
            break
        end
    end
    virus = virus(1:I);
    tts = tt(1:I);
    if I > 1
        virus = [virus(1),virus(I)];
        tts = [tt(1),tt(I)];
    end
    t = linspace(tts(1),tts(end));
    lb = [2.1e-7, 0];
    ub = [1,12];
    p0 = (lb+ub)./2;
    [vpfit resnorm] = lsqcurvefit(@paramfun1,p0,tts,virus,lb,ub);
%     array2table(vpfit)
%     %vpfit = [vpfit,]
%     %[vpfit resnorm] = lsqcurvefit(@paramfun1,vpfit,tts,virus,lb,ub);
     %[vpfit resnorm] = lsqcurvefit(@paramfun1,vpfit,tts,virus,lb,ub);
  %fit = paramfun1(vpfit,t);
 %   figure
   %   plot(t,fit,'r',tts,virus,'r*')
    % Virusfits(j,:) = vpfit; 
     %save('Virusfits.mat','Virusfits')
end