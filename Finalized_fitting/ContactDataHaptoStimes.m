function ret = ContactDataHaptoStimes(its,opt)
% function call for single instance of data fitting given a drawn infection
% start time

% inputs 
% its - which host 1-12 according to the order in the data tables
% opt - specify desired return
% opt = 1 -> parameter estimates and relative error
% opt = 2 -> generated dataset for Monte-Carlo simulations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import and organize data

% Innate response data (haptoglobin)
xs = readtable('HaptoDataContact.csv');

% Infection start times 
Starts = readtable('InfectionStartTimes.csv');
Starts = table2array(Starts);



% take log of Innate response to get everything on the same scale
xs = table2array(xs(:,3:end));
xs(xs>0) = log10(xs(xs>0));

% virus data
ses = readtable('ViremiaContact.csv');
ses = table2array(ses(:,3:end));

% adaptive data 
pvs = readtable('VNTContactData.csv');
pvs = table2array(pvs(:,3:end));


% days corresponding to data points (days since contact)
t2 = [0,2, 4, 6, 9,12,28];

t = t2;


% objective function 
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
function yy = paramfun1(p,t)
    

    lambda = (1/21.23).*min(x).*.8;
    k=p(1); d=(1/21.23); r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); % use disperse here ...
    x0 = [p(7),vv(1),p(8)]; % initial conditions 
    atilde=0;

     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(1:end,1)'.*sqwt3,yy(1:end,2)'.*sqwt,yy(:,3)'.*sqwt2];
end 




% function for linear regression on initial immune data
function ret = mdl(p,t)
    int = p(2); slope = p(1);
    ret = slope.*t + int;
end



% draw an infection start time 
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
        % redraw start time if too close to first positive measure 
        SIDX = randi([1 20000],1,1);
        Start = Starts(SIDX,j);
    end
    ind = find(t > Start);

    tTest = [Start,t(ind)];
    tt = tTest(1):.01:60;
    ttt = tTest(1):.001:30;

    time =   [t(ind(1)-1)  t(ind(1))];
    xtemp =   [min(x(ind(1)-1:ind(1)+1)) max(x(ind(1)-1:ind(1)))];
    pvtemp =   [min(pv(ind(1)-1:ind(1)+1)) max(pv(ind(1)-1:ind(1)))];


% fit regression models for initial conditions 
p0  = [abs((x(2)-x(1))./(t2(2)-t2(1))),1];
lb = [0, .01];
lu = [inf,inf];
[pfit_x resnorm] = lsqcurvefit(@mdl,p0,t2(1:3),x(1:3),lb,lu);
ff_x = @(t) pfit_x(1).*t + pfit_x(2);


p0  = [(pv(3)-pv(1))./(t2(3)-t2(1)),1];
lb = [0, .05];
lu = [inf,inf];
[pfit_pv resnorm] = lsqcurvefit(@mdl,p0,t2(1:3),pv(1:3),lb,lu);
ff_pv = @(t) pfit_pv(1).*t + pfit_pv(2);


% define new time series based on drawn infection start time and regression

    xNew = [ff_x(Start),xs(j,ind(1):end)];
    sNew = [0,ses(j,ind(1):end)];
    pvNew = [ff_pv(Start),pvs(j,ind(1):end)];
    s_idx = find(sNew > 0.5);
    s_idx = s_idx(1);
    
    % get estimates of inital viral load and viral growth rate
    vv = ViralGrowth(sNew,tTest-Start,s_idx);

    
    fprintf('-------------------------------------------------------------\n')

 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set non-linear least squares weights for each host, treated as hyper
% parameters 


sqwt = .1.*ones(1,length(pvNew));
s_ind = find(sNew > 0);
[MS IS] = max(sNew);
sqwt(s_ind(1):s_ind(end)+1) = 1;
if its ~= 10
sqwt(IS:IS+2) = 2;
elseif its == 10
sqwt(IS:IS+1) = 2;
end

sqwt2 = .1.*ones(1,length(pvNew));
pv_ind = find(pvNew > 0);
sqwt2(pv_ind(1):pv_ind(end)) = 2;
if its == 10
    sqwt2(end) = 5;
end
if its ~= 10
sqwt2(4) = 5;
sqwt2(5) = 3;
end

if its == 5
sqwt(4) = 0;
sqwt(3) = 2.25;
sqwt2=2;
sqwt3=1;


end

if its == 7
    if Start < t2(2)
    sqwt2(5) = 0;
    elseif Start >= t2(2)
        sqwt2(4) = 0;
    end
    
end
if its == 8
    if Start > t2(2)
    sqwt2(3) = .1;
    sqwt2(4) = 4;
    sqwt2(5) = 4;
    sqwt(3) = .1;
    sqwt(4) = 2;
    sqwt(1) = 2;
    elseif Start <= t2(2)
    sqwt2(4) = .1;
    sqwt2(5) = 4;
    sqwt2(6) = 4;
    sqwt(4) = .1;
    sqwt(5) = 2;
    sqwt(2) = 2;
    end
end

sqwt3 = .1.*ones(1,length(pvNew));
x_ind = find(xNew > 0);
sqwt3(x_ind(2):x_ind(end)) = 1;
sqwt3(end) = 1.5;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555





% set bound for regression

 lb = [.075,ceil(max(sNew)),0,0,min(2,vv(2)),.7*max(sNew),.8.*ff_x(Start),.8.*ff_pv(Start)];
 lu = [.5,11,2.25,1.5,15,1.3*max(sNew),1.2.*ff_x(Start),1.2.*ff_pv(Start)];


if its == 6 && Start > 2
    xNew(1) = x(2);

end

if vv(2) >= 6
    lb(5) = 2;
end

% set starting points for Newton-Rhapson 

p0 = (lb+lu)./2;
if its == 5
p0 = [0.2405    5.8776    2.5000    0.2813    3.9943    3.0000    3.9257    0.0859];
end

% fit data given drawn infection start times 
[pfit resnorm] = lsqcurvefit(@paramfun1,p0,tTest,[xNew(1:end).*sqwt3,sNew(1:end).*sqwt,pvNew(1:end).*sqwt2],lb,lu);

    
    
      lambda = (1/21.23).*min(x).*.8;
    k=pfit(1); d=(1/21.23); r=pfit(5); K=pfit(2); theta=0; delta=pfit(3); b=pfit(4); 
    x0 = [pfit(7),vv(1),pfit(8)]; % initial conditions 
    atilde=0; v0 = pfit(6);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [tplot,yy] = ode45(f,tt,x0);
   [tStore,fv] = ode45(f,tTest,x0);
   yyint = yy;
   tint = tplot;
   [mVi I] = max(yy(:,2));

  % get max and cum viral load 
Q = trapz(tint,yyint(:,2));
  cumviralContact(j) = Q;
    maxViralContact(j,:) = [mVi tt(I)];
    if opt == 1
     pfitsContactFinal = [pfit,vv(1),mVi,tt(I),Q,Start];
     ret = pfitsContactFinal;
    end



 
% set noise level (50%) and generate new datset for Monte-Carlo simulation 

xfitsC = fv(:,1)';
sfitsC = fv(:,2)';
sfitsC(sfitsC < 0) = -sfitsC(sfitsC < 0);
pvfitsC = fv(:,3)';
nLevel = [.5,.5,.5];

dd = GetData(xfitsC,sfitsC,pvfitsC,nLevel);
LL = length(xfitsC);
xNew2 = dd(1:LL);
sNew2 = dd(LL+1:2*LL);
pvNew2 = dd(2*LL+1:end);

% re estimate inital values for A, I

p0  = [abs((xNew2(2)-xNew2(1))./(t2(2)-t2(1))),1];
lb = [0, .01];
lu = [inf,inf];
[pfit_x resnorm] = lsqcurvefit(@mdl,pfit_x,t2(1:3),GetDataInitial(ff_x(t2(1:3)),.5),lb,lu);
ff_x = @(t) pfit_x(1).*t + pfit_x(2);


p0  = [(pvNew2(3)-pvNew2(1))./(t2(3)-t2(1)),1];
lb = [0, .05];
lu = [inf,inf];
[pfit_pv resnorm] = lsqcurvefit(@mdl,pfit_pv,t2(1:3),GetDataInitial(ff_pv(t2(1:3)),.5),lb,lu);
ff_pv = @(t) pfit_pv(1).*t + pfit_pv(2);

lu(2) = min(9,1.25*max([sNew,sNew2]));

% get P0 and inital estimate of viral growth rate for generated dataset
vvNew = ViralGrowth(sNew2,tTest-Start,s_idx);
 lb(end-1) = 1e-4;
    lu(end-1) = 2*min(max(s),max(sNew2));
vvOld = vv;
vv = vvNew;
lb(end) = .8*min(pvNew(1),pvNew2(1));  
lu(end) = 1.2*max(pvNew(1),pvNew2(1));  

% fit the generated dataset 
 lb = [.075,ceil(min([max(sNew),max(sNew2)])),0,0,min([2,vvOld(2),vv(2)]),.7.*min(max(sNew),max(sNew2)),min(.8.*ff_x(Start),.8.*pfit(7)),min(.8.*ff_pv(Start),.8.*pfit(8))];
 lu = [.5,11,2.25,1.5,15,1.3*max(max(sNew),max(sNew2)),max(1.2.*ff_x(Start),1.2.*pfit(7)),max(1.2.*ff_pv(Start),1.2.*pfit(8))];

 [pfit2 resnorm] = lsqcurvefit(@paramfun1,pfit,tTest,[xNew2(1:end).*sqwt3,sNew2(1:end).*sqwt,pvNew2(1:end).*sqwt2],lb,lu);
 lambda = (1/21.23).*min(x).*.8;
    k=pfit2(1); d=(1/21.23); r=pfit2(5); K=pfit2(2); theta=0; delta=pfit2(3); b=pfit2(4); % use disperse here ...
    x0 = [pfit2(7),vv(1),pfit2(8)]; % initial conditions 
    atilde=0; v0 = pfit2(6);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
  % recalculate cum and max viral 
   [tplot,yy] = ode45(f,tt,x0);its
   [tStore,fv] = ode45(f,tTest,x0);
   yyint = yy;
   tint = tplot;
   [mVi I] = max(yy(:,2));

Q = trapz(tint,yyint(:,2));
pfitsContactFinal2 = [pfit2,vv(1),mVi,tt(I),Q,Start];

if opt == 1
RE = abs(pfitsContactFinal2-pfitsContactFinal)./pfitsContactFinal;
ret = [pfitsContactFinal;pfitsContactFinal2;RE];
end
if opt == 2
    ret = [xNew2;sNew2;pvNew2;tTest];
end
end


end