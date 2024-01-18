function [paramfit] = FittingforCIHapto(dd,index)
% input is synthetic dataset and index for appropriate dataset baseline
% parameter values
points = length(dd)./3;
x = dd(1:points);
s = dd(points+1:2*points);
vv = ViralGrowth(s);
pv = dd(2*points+1:end);

params = readtable('FitParamsHapto.csv');
params = table2array(params(:,4:10));
p0 = params(index,:);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
tt = t(1):.01:30;
ttt = t(1):.001:30;
% parameter boundries 

 sqwt = 1; %r (virus)
sqwt2 = 1; % g (adapt)
sqwt3 = 1; % b (innate)
%lb = [0,4,0,0,0,vv(2),vv(1)-5*eps,9e-7];
%lu = [.5,10,2.0,1,0,5,max(s),9e-3];
lb = [.04,3.75,0,0,0,vv(2),.9*vv(1)];
lu = [.5,11,3,1,0,6,1.1*vv(1)];
%array2table([lb;lu])
lb(end-2) = .8*x(1);
lu(end-2) = 1.0*x(1);
if lu(end-2) < 1
    lb(end-2) = 1.5;
    lu(end-2) = 2.0;
end

options=odeset('RelTol',1e-12,'AbsTol',1e-12);
% objective function 
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
   yy = [yy(1:end,1)'*sqwt3,yy(1:end,2)'*sqwt,yy(:,3)'*sqwt2];
end 


   [pfit resnorm] = lsqcurvefit(@paramfun1,p0,t,[x*sqwt3,s*sqwt,pv*sqwt2],lb,lu);
   %[pfit resnorm] = lsqcurvefit(@paramfun1,pfit,t,[x*sqwt3,s*sqwt,pv*sqwt2],lb,lu);
         lambda1=(1/11.8)*pfit(end-3);
    
    
        lambda2 = (1/21.23)*x(end);
        lambda = min(lambda1,lambda2);
    k=pfit(1); d=1/21.23; r=pfit(6); K=pfit(2); theta=0; delta=pfit(3); b=pfit(4); % use disperse here ...
    x0 = [pfit(5),vv(1),1e-4*pfit(5)]; % initial conditions 
    atilde=0; v0 = pfit(7);
  f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [tplot,yy] = ode45(f,tt,x0);
   [tint,yyint] = ode45(f,ttt,x0);
   [mVi I] = max(yy(:,2));

   Q = trapz(tint,yyint(:,2));
  

    paramfit = [pfit,vv(1),mVi,tt(I),Q];
   %paramfit = [vv(2),vv(1),pfit];
end



