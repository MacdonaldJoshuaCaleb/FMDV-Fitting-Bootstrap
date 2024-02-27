function vpfit = ViralGrowth(ses,tt,s_idx)
% input 
% viral load data - ses
% observation times - tt
% index of first positive viral load measure - s_idx

% return estiamte for P0 and r

% objective function, assumes simple exponential growth
function yy = paramfun1(p,t)
    x0 = p(1);
    r = p(2);
    
    f = @(t) x0.*exp(r.*t);
   yy = f(t);
end 

% truncate to only use up to first postive measured load 
virus = ses;
    I = s_idx;
    virus = virus(1:I);
    tts = tt(1:I);
    if I > 1
        virus = [virus(1),virus(I)];
        tts = [tt(1),tt(I)];
    end
    % define time stepping, lower, and upper bounds, and initial guess
    t = linspace(tts(1),tts(end));
    lb = [.005, 1];
    ub = [1,12];
    p0 = [.01,4];
    
    % fit the data 
    [vpfit resnorm] = lsqcurvefit(@paramfun1,p0,tts,virus,lb,ub);


end