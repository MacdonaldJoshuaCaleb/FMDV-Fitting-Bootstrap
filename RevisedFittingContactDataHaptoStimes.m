function ret = RevisedFittingContactDataHaptoStimes(its,opt)
%close all
% turn figure display on  or off
set(0,'DefaultFigureVisible','on')
% import and organize data
xs = readtable('HaptoDataContact.csv');

Starts = readtable('InfectionStartTimes.csv');
Starts = table2array(Starts);
%Starts = mean(Starts);

groups = xs(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xs(:,2);
IDs = table2array(IDs);
%IDs = string(groups);

xs = table2array(xs(:,3:end));
xs(xs>0) = log10(xs(xs>0));

ses = readtable('ViremiaContact.csv');
ses = table2array(ses(:,3:end));

pvs = readtable('VNTContactData.csv');
pvs = table2array(pvs(:,3:end));

% set up solution bin to store fit parameters for analysis 


% initial guesss
% initial conditions 
%x0 = [2.4223;0.8164;2];
% params 
%lambda=5;k=1;d=0.5;r=8;K=8;theta=0.5;delta=0.7;atilde=0.1;b=0.0001;
params = readtable('FitParamsHapto.csv');
params = table2array(params(:,4:10));

% k K delta b I0 r PS P0 
% days corresponding to data points 
t2 = [0,2, 4, 6, 9,12,28];
% t = [0,2, 4, 6, 8, 11, 14, 30];
t = t2;

% tt = t(1):.01:30;
% ttt = t(1):.001:30;
% parameter boundries 

% sqwt = 3; %r (virus)
% sqwt2 = 3; % g (adapt)
% sqwt3 = 1; % b (innate)
% objective function 
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
function yy = paramfun1(p,t)
    
    lambda1=(1/21.23)*p(end-2);
    
    
        lambda2 = (1/21.23)*x(end);
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/21.23; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),vv(1),1e-4*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(1:end,1)'.*sqwt3,yy(1:end,2)'.*sqwt,yy(:,3)'.*sqwt2];
end 





%its = 1:12;
for jj = its:its
    j = its;
    x = xs(j,:);
    s = ses(j,:);
    pv = pvs(j,:);
    SIDX = randi([1 20000],1,1);
    Start = Starts(SIDX,j);
    for kk = 1:length(s)
        if s(kk) > 0
            ind2 = t(kk);
            break
        end
    end
    while Start >= ind2-1
        SIDX = randi([1 20000],1,1);
        Start = Starts(SIDX,j);
    end
    ind = find(t > Start);

    tTest = [Start,t(ind)];
    tt = tTest(1):.01:30;
    ttt = tTest(1):.001:30;

    time =   [t(ind(1)-1)  t(ind(1))];
    xtemp =   [min(x(ind(1)-1:ind(1)+1)) max(x(ind(1)-1:ind(1)))];
%     stemp = [s(ind(1)-1) s(ind(1))];
%     pvtemp = [pv(ind(1)-1) pv(ind(1))];
    
    xStart = interp1(time,xtemp,Start,'linear');
%     
%     if xtemp(end) < max(x(1:3))
%         xStart = xStart = interp1(time,xtemp,Start,'linear');
%     end
% %     sStart = interp1(time,stemp,Start,'linear');
%    pvStart = interp1(time,pvtemp,Start,'linear');
    xNew = [xStart,xs(j,ind(1):end)];
    sNew = [0,ses(j,ind(1):end)];
    pvNew = [0,pvs(j,ind(1):end)];

    vv = ViralGrowth(sNew,tTest);
    
    fprintf('-------------------------------------------------------------\n')
% % % % % if its==8
% % % % %  sqwt = 1.*ones(1,length(sNew)); %r (virus)
% % % % %  [M ind] = max(sNew);
% % % % %  sqwt(ind) = 8;
% % % % %  sqwt(ind+1:ind+2) = 2;
% % % % %  pvi = find(pvNew > 0);
% % % % %    sqwt2 = 1.*ones(1,length(pvNew));
% % % % %    sqwt(pvi(1)) = .25;
% % % % %     sqwt3 = .5;
% % % % %  
% % % % % end
%  sqwt(ind(1)) = 8;
%  sqwt(ind(end)) = 8;
% sqwt2 = 1.*ones(1,length(pvNew)); % g (adapt)
% ind2 = find(pvNew > 0);
% sqwt2(end-2:end) = 2;
% sqwt3 = 1.*ones(1,length(xNew)); % b (innate)
% [M I] = max(xNew);
% sqwt3(end-1:end) = 2;
%if its == 4
%     sqwt = 1.*ones(1,length(sNew));
%     ind = find(sNew > 0);
%     sqwt(ind(1)-1:ind(end)+1) = 2;
%     %sqwt(ind(1)) = 1;
%     sqwt2 = 1*ones(1,length(pvNew));
%     sqwt2(end-1:end) = 2;
%     sqwt3 = 1;
%end
% if its == 7
%     sqwt = 1.*ones(1,length(sNew));
%     [M ind] = max(sNew);
%     sqwt(ind:ind+1) = 2;
%     sqwt(ind-1) = 4;
%     sqwt2 = 1;
%     sqwt3 = 1;
% end
%if its == 9
% % % % % if its ~= 8
    sqwt = 1.*ones(1,length(sNew));
    [M ind] = max(sNew);
    sqwt(ind:ind+1) = 5;
    sqwt(ind-1) = 10;
    sqwt2 = 1;
    sqwt3 = 1.*ones(1,length(xNew));
    [M xI] = max(xNew);
    sqwt3(xI) = 3;
%     if its == 8
%         pvi = find(pvNew > 0);
%         sqwt2 = 1.*ones(1,length(pvNew));
%         sqwt2(pvi(1:end-1)) = 0;
%         sqwt2(pvi(1)-1) = 2;
%         sqwt = 1.*ones(1,length(sNew));
%         sqwt(ind) = 10;
%         sqwt(ind+1) = .75;
%         %sqwt(ind-1) = 4;
%     end
    
    if its ~= 7
        pvi = find(pvNew > 0);
        sqwt2 = 1.*ones(1,length(pvNew));
        sqwt2(pvi(1:end-1)) = 0;
        if its == 8
            sqwt2(pvi(2:end)) = 2;
        end
        sqwt2(pvi(1)-1) = 2;
        sqwt = 1.*ones(1,length(sNew));
        sqwt(ind) = 10;
        sqwt(ind+1) = .5;
        %sqwt(ind-1) = 4;
    end
%end
if its == 11
        sqwt = 1.*ones(1,length(sNew));
    [M ind] = max(sNew);
    sqwt(ind:ind+1) = 5;
    sqwt(ind-1) = 10;
end
if its == 9
     sqwt = 1.*ones(1,length(sNew));
    [M ind] = max(sNew);
    sqwt(ind:ind+1) = 5;
    sqwt(ind-1) = 5;
    pvi = find(pvNew > 0);
        sqwt2 = 1.*ones(1,length(pvNew));
        sqwt2(pvi(1)) = 3;
        sqwt2(pvi(1)-1) = 3;
        
end
% 
% if its == 11
%     sqwt = 1.*ones(1,length(sNew));
%     [M ind] = max(sNew);
%     sqwt(ind:ind+1) = 1;
%     sqwt(ind-1) = 1;
%     sqwt2 = 1;
%     sqwt3 = 1;
% end
[SM xind] = max(sNew);
sNew(xind:end);
sind = find(sNew > 0);
% lb = [0.08,0,0,0,0,vv(2),.3*max(s)];
% lu = [.6,min(1.25*max(sNew),9),4.5,1.5,0,11,.55*max(s)];


 lb = [0.04,0,0,0,0,vv(2),.8*vv(1)];
 lu = [.5,min(1.25*max(sNew),9),4.75,1.5,0,11,1.1*vv(1)];

if its >= 8
    lb(end) = .9*vv(1);
    lu(end) = 1.1.*vv(1);
end
if its == 6 && Start > 2
    xNew(1) = x(2);
    %lb(end) = .9*vv(1);
    %lu(end) = 1.1.*vv(1);
end

if vv(2) >= 7
    lb(6) = 3;
end
lb(end-2) = .8*xNew(1);
%lu(end-2) = min([xNew(1:4),3.5]);
lu(end-2) = 1.2*xNew(1);
% sqwt = 1;
% sqwt2 = 1;
% sqwt3 = 1;
p0 = mean([lb;lu]);
%p0(end) = 5;
%p0(3) = 2;
p0(end-2) = .9*min(xNew(1:3));
p0(end-1) = vv(2);
p0(1) = .1;
p0(3) = .4;

   [pfit resnorm] = lsqcurvefit(@paramfun1,p0,tTest,[xNew(1:end).*sqwt3,sNew(1:end).*sqwt,pvNew(1:end).*sqwt2],lb,lu);
   %[pfit resnorm] = lsqcurvefit(@paramfun1,pfit,tTest,[xNew(1:end).*sqwt3,sNew(1:end).*sqwt,pvNew(1:end).*sqwt2],lb,lu);
% array2table(pfit)
%    bFits = paramfun1(pfit,t);
%    xfitsC = bFits(1:length(t))./sqwt3;
%    sfitsC = bFits(length(t)+1:2*(length(t)))./sqwt;
%    pvfitsC = bFits(2*(length(t))+1:end)./sqwt2;
%    %pfitsContactFinal(j,:) = pfit;   %lambda1=(1/11.8)*pfit(end-2);
    
    
        lambda1=(1/21.23)*pfit(end-2);
    
    
        lambda2 = (1/21.23)*x(end);
        lambda = min(lambda1,lambda2);
    k=pfit(1); d=1/21.23; r=pfit(6); K=pfit(2); theta=0; delta=pfit(3); b=pfit(4); % use disperse here ...
    x0 = [pfit(5),vv(1),1e-4*pfit(5)]; % initial conditions 
    atilde=0; v0 = pfit(7);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [tplot,yy] = ode45(f,tt,x0);
   [tStore,fv] = ode45(f,tTest,x0);
   yyint = yy;
   tint = tplot;
   [mVi I] = max(yy(:,2));
%    fits = paramfun1(pfit,tt);
%    xfit = fits(1:length(tt));
%    sfit = fits(length(tt)+1:2*(length(tt)));
%    [mVi I] = max(sfit/sqwt);
Q = trapz(tint,yyint(:,2));
  cumviralContact(j) = Q;
    maxViralContact(j,:) = [mVi tt(I)];
    if opt == 1
     pfitsContactFinal = [pfit,vv(1),mVi,tt(I),Q,Start];
     ret = pfitsContactFinal;
    end
%    pvfit = fits(2*length(tt)+1:end);
    str = strcat('ID',' ',num2str(IDs(j)),', ',extractBetween(groups(j),1,4));
    str2= strcat('ID',num2str(IDs(j)),groups(j),'A_0','Hapto');
%    length(tt(1:end-1))
%    length(xfit)
   % figure
% % % % %    hold on
% % % % %    plot(tt(1:end),yy(:,1),'b',tt(1:end),yy(:,2),'r',tt,yy(:,3),'g','linewidth',2)
% % % % %    plot(t2,x,'b*',t2,s,'r+',t2,pv,'go','MarkerSize',20,'linewidth',2)
% % % % %    xline(Start,'k--','linewidth',2)
% % % % %    hold off
% % % % %    ylim([0 12])
% % % % %    xlim([0 15])
% % % % %   if j == 1
% % % % %   % legend('I(\tau)','P(\tau)','A(\tau)','log_{10}(I_{SAA})','FMDV','A_{VNT}')
% % % % %    end
% % % % %    title(str)
% % % % %    set(gca,'FontSize',16)
%     baseFileName = sprintf(str2);
%         fname = '~/Documents/MATLAB/Model_identifiability/Plots';
%      saveas(gca, fullfile(fname, baseFileName), 'png');
% % % % % %    dvdtau = (gradient(sfit(:)/sqwt)./gradient(tt(1:end-1)));
 %   [mindV, index] = min(dvdtau);
  %  for w = 1:length(sfit)
   %     if abs(dvdtau(w)) <= .05
    %        if w > index
     %       VTC(j) = tt(w);
      %      break
 %           end
  %      end
  %  end
  
xfitsC = fv(:,1)';
sfitsC = fv(:,2)';
sfitsC(sfitsC < 0) = -sfitsC(sfitsC < 0);
pvfitsC = fv(:,3)';
nLevel = [.4,.4,.4];
dd = GetData(xfitsC,sfitsC,pvfitsC,nLevel);
LL = length(xfitsC);
xNew2 = dd(1:LL);
sNew2 = dd(LL+1:2*LL);
pvNew2 = dd(2*LL+1:end);
lu(end-2) = 1.2*max(xNew2(1),xNew(1));
lb(end-2) = .8*min(xNew2(1),xNew(1));
% lu(end) = .55*max(max(sNew2),max(sNew));
% lb(end) = .3*min(max(sNew2),max(sNew));

lu(2) = min(9,1.25*max([sNew,sNew2]));
vvNew = ViralGrowth(sNew2,tTest);
%lu(end-1) = 2*max(vv(2),vvNew(2));
 lb(end) = .8*max(vv(1),vvNew(1));
    lu(end) = 1.2*max(vv(1),vvNew(1));

vv = vvNew;
   


%lb(end-1) = min(vv(2),pfit(6));
if lb(end-1) >= 7
    lb(6) = 3;
end

%         array2table([lb;lu])
%     if its ~= 7
%         pvi = find(pvNew > 0);
%         sqwt2 = 1.*ones(1,length(pvNew));
%         sqwt2(pvi(1:end-1)) = 0;
%         sqwt2(pvi(1)-1) = 2;
%         sqwt = 1.*ones(1,length(sNew));
%         sqwt(ind) = 10;
%         sqwt(ind+1) = .75;
%         %sqwt(ind-1) = 4;
%     end

[pfit2 resnorm] = lsqcurvefit(@paramfun1,pfit,tTest,[xNew2(1:end).*sqwt3,sNew2(1:end).*sqwt,pvNew2(1:end).*sqwt2],lb,lu);
    lambda1=(1/21.23)*pfit2(end-2);
    
    
        lambda2 = (1/21.23)*x(end);
        lambda = min(lambda1,lambda2);
    k=pfit2(1); d=1/21.23; r=pfit2(6); K=pfit2(2); theta=0; delta=pfit2(3); b=pfit2(4); % use disperse here ...
    x0 = [pfit2(5),vv(1),1e-4*pfit2(5)]; % initial conditions 
    atilde=0; v0 = pfit2(7);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [tplot,yy] = ode45(f,tt,x0);
   [tStore,fv] = ode45(f,tTest,x0);
   yyint = yy;
   tint = tplot;
   [mVi I] = max(yy(:,2));
%    fits = paramfun1(pfit,tt);
%    xfit = fits(1:length(tt));
%    sfit = fits(length(tt)+1:2*(length(tt)));
%    [mVi I] = max(sfit/sqwt);
Q = trapz(tint,yyint(:,2));
pfitsContactFinal2 = [pfit2,vv(1),mVi,tt(I),Q,Start];
if opt == 1
RE = abs(pfitsContactFinal2-pfitsContactFinal)./pfitsContactFinal;
ret = [pfitsContactFinal;pfitsContactFinal2;RE];
end
if opt == 2
    ret = [xfitsC;sfitsC;pvfitsC];
end
end
% array2table(pfitsContactFinal)
% save('pfitsContactFinalHapto.mat','pfitsContactFinal')
% save('maxViralContactFinalHapto.mat','maxViralContact')
% save('cumviralContactFinalHapto.mat','cumviralContact')
%save('VTCFinal.mat','VTC')
% save('xfitsCFinalHapto.mat','xfitsC')
% save('sfitsCFinalHapto.mat','sfitsC')
% save('pvfitsCFinalHapto.mat','pvfitsC')

end