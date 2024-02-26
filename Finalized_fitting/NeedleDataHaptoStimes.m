function ret = NeedleDataHaptoStimes(its,opt)
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
% import and organize data

% Haptoglobin data 
xs = readtable('HaptoDataNeedle.csv');

% grab out info for plot title generation 
groups = xs(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xs(:,2);
IDs = table2array(IDs);

% log_10 transform the innate data to be on same sacle as others
xs = table2array(xs(:,3:end));
xs(xs>0) = log10(xs(xs>0));

% read in virus data
ses = readtable('ViremiaNeedle.csv');
ses = table2array(ses(:,3:end));

% read in VNT data
pvs = readtable('VNTNeedleData.csv');
pvs = table2array(pvs(:,3:end));

% time in experiment days (since needle infected)
t = [0,2, 4, 6, 8, 11, 14, 30];
% store a copy for later calculations 
t2 = t;

% objective function 
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
function yy = paramfun1(p,t)
    
 
       lambda = (1/21.23).*min(x).*.8; % fix lambda for each host so that equib. is 80% minimum data value
       % define paramters to be inferred 
    k=p(1); d=(1/21.23); r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); 
    % vv refers to initial estiamte of initial viral load form viralgrowth .m later in the code 
    x0 = [p(7),vv(1),p(8)]; % initial conditions 
    atilde=0;
    % the ode
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   % return concentatnated array of weighted model solutions 
   yy = [yy(1:end,1)'.*sqwt3,yy(1:end,2)'.*sqwt,yy(:,3)'.*sqwt2];
end 




% regression for initial innate/adaptive values as defined 
function ret = mdl(p,t)
    int = p(2); slope = p(1);
    ret = slope.*t + int;
end

function ret = mdl2(p,t)
    ret = p(1).*exp(p(2).*t);
end

for jj = its:its 
    j = its;
    x = xs(j,:);
    s = ses(j,:);
    pv = pvs(j,:);
    SIDX = randi([1 20000],1,1);
    Start = 0;
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
    tt = tTest(1):.01:60;
    ttt = tTest(1):.001:30;

    time =   [t(ind(1)-1)  t(ind(1))];
    xtemp =   [min(x(ind(1)-1:ind(1)+1)) max(x(ind(1)-1:ind(1)))];
    pvtemp =   [min(pv(ind(1)-1:ind(1)+1)) max(pv(ind(1)-1:ind(1)))];
%     stemp = [s(ind(1)-1) s(ind(1))];
%     pvtemp = [pv(ind(1)-1) pv(ind(1))];


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



    xNew = [ff_x(Start),xs(j,ind(1):end)];
    sNew = [0,ses(j,ind(1):end)];
    pvNew = [ff_pv(Start),pvs(j,ind(1):end)];
    s_idx = find(sNew > 0.5);
    s_idx = s_idx(1);
    vv = ViralGrowth(sNew,tTest-Start,s_idx);

    
    fprintf('-------------------------------------------------------------\n')

 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_s = length(find(s>0));
len_x = length(find(x>0));
len_pv = length(find(pv>0));
sqwt = .1.*ones(1,length(pvNew));
s_ind = find(sNew > 0);
[MS IS] = max(sNew);
sqwt(s_ind(1):s_ind(end)+1) = 1;

sqwt(IS:IS+2) = 2;


sqwt2 = .1.*ones(1,length(pvNew));
pv_ind = find(pvNew > 0);
sqwt2(pv_ind(1):pv_ind(end)) = 2;


sqwt2(4) = 5;
sqwt2(5) = 3;





sqwt3 = 1*ones(1,length(pvNew));
x_ind = find(xNew > 0);
sqwt3(x_ind(2):x_ind(end)) = 1;
sqwt3(end) = 0;
% sqwt3(end-1) = 
sqwt3(2) = 3;


% sqwt = 1;
% sqwt2 = 1;
% sqwt3 = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555


% [MM II] = max(xNew);
% sqwt3(II) = 3;
% sqwt3(end-3:end)= 3;
[SM xind] = max(sNew);
sNew(xind:end);
sind = find(sNew > 0);
% lb = [0.08,0,0,0,0,vv(2),.3*max(s)];
% lu = [.6,min(1.25*max(sNew),9),4.5,1.5,0,11,.55*max(s)];


 lb = [.075,ceil(max(sNew)),0,0,min(2,vv(2)),.7*max(sNew),min(.8*x(1),2),.8.*ff_pv(Start)];
 lu = [1,11,2.25,1.5,15,1.3*max(sNew),max(2,1.2*x(1)),1.2.*ff_pv(Start)];

% if its >= 8
%     lb(end) = .9*vv(1);
%     lu(end) = 1.1.*vv(1);
% end
if its == 6 && Start > 2
    xNew(1) = x(2);
    %lb(end) = .9*vv(1);
    %lu(end) = 1.1.*vv(1);
end

if vv(2) >= 6
    lb(5) = 2;
end


% p0(end) = 1/20;
% params = readtable('params_may3.csv');
% params = table2array(params);
% p0 = params(its,1:8);
p0 = (lb+lu)./2;
if its == 5
p0 = [0.2405    5.8776    2.5000    0.2813    3.9943    3.0000    3.9257    0.0859];
end

% p0(end-1) = 5;
% if p0(6) > 5
%     p0(6) = 2;
% end

   [pfit resnorm] = lsqcurvefit(@paramfun1,p0,tTest,[x(1:end).*sqwt3,s(1:end).*sqwt,pv(1:end).*sqwt2],lb,lu);

    
    
      lambda = (1/21.23).*min(x).*.8;
    k=pfit(1); d=(1/21.23); r=pfit(5); K=pfit(2); theta=0; delta=pfit(3); b=pfit(4); % use disperse here ...
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

  
xfitsC = fv(:,1)';
sfitsC = fv(:,2)';
sfitsC(sfitsC < 0) = -sfitsC(sfitsC < 0);
pvfitsC = fv(:,3)';
nLevel = [.5,.5,.5];
% figure
% hold on
% for zz = 1:200
dd = GetData(xfitsC,sfitsC,pvfitsC,nLevel);
LL = length(xfitsC);
xNew2 = dd(1:LL);
sNew2 = dd(LL+1:2*LL);
pvNew2 = dd(2*LL+1:end);


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

vvNew = ViralGrowth(sNew2,tTest-Start,s_idx);
%lu(end-1) = 2*max(vv(2),vvNew(2));
 lb(end-1) = 1e-4;
    lu(end-1) = 2*min(max(s),max(sNew2));
vvOld = vv;
vv = vvNew;
lb(end) = .8*min(pvNew(1),pvNew2(1));  
lu(end) = 1.2*max(pvNew(1),pvNew2(1));  







 lb = [.075,ceil(min([max(sNew),max(sNew2)])),0,0,min([2,vvOld(2),vv(2)]),.7.*min(max(sNew),max(sNew2)),min([.8.*pfit(7),2,.8*xNew2(1)]),min(.8.*ff_pv(Start),.8.*pfit(8))];
 lu = [1,11,2.25,1.5,15,1.3*max(max(sNew),max(sNew2)),max(1.2*x(1),1.2*xNew2(1)),max(1.2.*ff_pv(Start),1.2.*pfit(8))];

 
 
[pfit2 resnorm] = lsqcurvefit(@paramfun1,pfit,tTest,[xNew2(1:end).*sqwt3,sNew2(1:end).*sqwt,pvNew2(1:end).*sqwt2],lb,lu);
 lambda = (1/21.23).*min(x).*.8;
    k=pfit2(1); d=(1/21.23); r=pfit2(5); K=pfit2(2); theta=0; delta=pfit2(3); b=pfit2(4); % use disperse here ...
    x0 = [pfit2(7),vv(1),pfit2(8)]; % initial conditions 
    atilde=0; v0 = pfit2(6);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [tplot,yy] = ode45(f,tt,x0);
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
