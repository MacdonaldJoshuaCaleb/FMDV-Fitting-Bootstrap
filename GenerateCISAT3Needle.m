function [] = GenerateCISAT3Needle
close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsN = readtable('pvfitsN.csv');
pvfitsN = table2array(pvfitsN);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsN = readtable('sfitsN.csv');
sfitsN = table2array(sfitsN);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsN = readtable('xfitsN.csv');
xfitsN = table2array(xfitsN);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsN = readtable('SAADataNeedle.csv');
xsN = table2array(xsN(:,3:end));
xsN(xsN>0) = log10(xsN(xsN>0));
sesN = readtable('ViremiaNeedle.csv');
sesN = table2array(sesN(:,3:end));
pvsN = readtable('VNTNeedleData.csv');
pvsN = table2array(pvsN(:,3:end));
% xsN = readtable('SAADataNeedle.csv');
% xsN = table2array(xsN(:,3:end));
% xsN(xsN>0) = log10(xsN(xsN>0));
% sesN = readtable('ViremiaNeedle.csv');
% sesN = table2array(sesN(:,3:end));
% pvsN = readtable('VNTNeedleData.csv');
% pvsN = table2array(pvsN(:,3:end));
params = readtable('FitParamsNeedle.csv');
params = table2array(params(:,4:end));
tt = 0:.1:30;
its = 10000;

Fits9 = zeros(its,3*length(tt));
 Fits10 = zeros(its,3*length(tt));
Fits11 = zeros(its,3*length(tt));
Fits12 = zeros(its,3*length(tt));

params9 = zeros(its,11);
 params10 = zeros(its,11);
params11 = zeros(its,11);
params12 = zeros(its,11);

%cumviralContact1 = zeros(its,1);
% cumviralContact2 = zeros(its,1);
% cumviralContact3 = zeros(its,1);
% cumviralContact4 = zeros(its,1);
% 
% maxViralContact1 = zeros(its,2);
% maxViralContact2 = zeros(its,2);
% maxViralContact3 = zeros(its,2);
% maxViralContact4 = zeros(its,2);

SData9 = zeros(its,24);
 SData10 = zeros(its,24);
SData11 = zeros(its,24);
SData12 = zeros(its,24);

function yy = paramfun1(p,t)
    vv = ViralGrowth(dd(length(xfitsN(index,:))+1:2*length(xfitsN(index,:))));
    lambda1=(1/11.18)*p(end-3);
    
    
        lambda2 = (1/11.18)*dd(length(xfitsN(index,:)));
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/11.18; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),vv(1),p(8)*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(1:end,1)',yy(1:end,2)',yy(:,3)'];
end 
in = [9,10,11,12];
for count = 1:its
for jj = 1:4
    index = in(jj);
nLevel = [.4,.4,.4];
dd = GetData(xfitsN(index,:),sfitsN(index,:),pvfitsN(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 9
    params9(count,:) = FittingforCINeedle(dd,index);
    Fits9(count,:) = paramfun1(params9(count,:),tt);
    %[mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    %cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    %maxViralContact1(count,:) = [mVi tt(I)];
    SData9(count,:) = dd;
end
if index == 10
    params10(count,:) = FittingforCINeedle(dd,index);
    Fits10(count,:) = paramfun1(params10(count,:),tt);
    %[mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    %cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    %maxViralContact2(count,:) = [mVi tt(I)];
    SData10(count,:) = dd;
end
if index == 11
    params11(count,:) = FittingforCINeedle(dd,index);
    Fits11(count,:) = paramfun1(params11(count,:),tt);
    %[mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    %cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    %maxViralContact3(count,:) = [mVi tt(I)];
    SData11(count,:) = dd;
end
if index == 12
    params12(count,:) = FittingforCINeedle(dd,index);
    Fits12(count,:) = paramfun1(params12(count,:),tt);
    %[mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    %cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    %maxViralContact4(count,:) = [mVi tt(I)];
    SData12(count,:) = dd;
end

end
SAT3means(count,:) = (params9(count,:)+params10(count,:)+params11(count,:)+params12(count,:))./4;
SAT3devs(count,:) = std([params9(count,:);params10(count,:);  params11(count,:); params12(count,:)]);
RE9(count,:) = abs(params9(count,:)-params(9,:))./params(9,:);
RE10(count,:) = abs(params10(count,:)-params(10,:))./params(10,:);
RE11(count,:) = abs(params11(count,:)-params(11,:))./params(11,:);
RE12(count,:) = abs(params12(count,:) - params(12,:))./params(12,:); 
% f = @(a,index) 4.*(log(a)-log(SAT1means(count,index))-psi(a) + mean(log([params1(count,index); params2(count,index); params3(count,index); params4(count,index)])));
% a0 = .5;
% for ss = 1:11
% a = a0 - f(a0,ss)./((f(a0+.01,ss)-f(a0,ss))./.01);
% temp = a0;
%     while abs(a - temp) > 1e-5
%         temp = a;
%         a = temp - f(temp,ss)./((f(temp+.01,ss)-f(temp,ss))./.01);
%         if a < 0
%             a = temp;
%             break
%         end
%     end
% alphas(ss) = a;
% betas(ss) = a./SAT1means(count,ss);
% alphasSAT1(count,ss) = alphas(ss);
% betasSAT1(count,ss) = betas(ss);
%end

fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 3*count,3*its)
fprintf('-----------------------\n')
end


save('Fits9Needle.mat','Fits9')
 save('Fits10Needle.mat','Fits10')
save('Fits11Needle.mat','Fits11')
save('Fits12Needle.mat','Fits12')

save('params9Needle.mat','params9')
 save('params10Needle.mat','params10')
save('params11Needle.mat','params11')
save('params12Needle.mat','params12')

ARE9 = mean(RE9)*100;
ARE10 = mean(RE10)*100;
ARE11 = mean(RE11)*100;
ARE12 = mean(RE12)*100;
array2table(ARE9)
array2table(ARE10)
array2table(ARE11)
array2table(ARE12)
for ss = 1:11
[pHat,pCI] = lognfit(params9(:,ss));
pHats9(ss,:) = pHat;
pCImu9(ss,:) = pCI(:,1);
pCIsig9(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params10(:,ss));
pHats10(ss,:) = pHat;
pCImu10(ss,:) = pCI(:,1);
pCIsig10(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params11(:,ss));
pHats11(ss,:) = pHat;
pCImu11(ss,:) = pCI(:,1);
pCIsig11(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params12(:,ss));
pHats12(ss,:) = pHat;
pCImu12(ss,:) = pCI(:,1);
pCIsig12(ss,:) = pCI(:,2);
end

%array2table(alphas)
%array2table(betas)
for ss = 1:11
    figure
    hold on
    %histogram(params1(:,ss))
    x = .4*min(params9(:,ss)):.0001:1.6*max(params9(:,ss));
    plot(x,lognpdf(x,pHats9(ss,1),pHats9(ss,2)),'linewidth',2)
    hold off
    xlim([min(x) max(x)])
end
% array2table(alphasSAT1)
% array2table(betasSAT1)
save('SData9Needle.mat','SData9')
save('SData10Needle.mat','SData10')
save('SData11Needle.mat','SData11')
save('SData12Needle.mat','SData12')
save('SAT3meansNeedle.mat','SAT3means')
save('SAT3devsNeedle.mat','SAT3devs')
save('pHats9Needle.mat','pHats9')
save('pCImu9Needle.mat','pCImu9')
save('pCIsig9Needle.mat','pCIsig9')
 save('pHats10Needle.mat','pHats10')
save('pCImu10Needle.mat','pCImu10')
save('pCIsig10Needle.mat','pCIsig10')
save('pHats11Needle.mat','pHats11')
save('pCImu11Needle.mat','pCImu11')
save('pCIsig11Needle.mat','pCIsig11')
save('pHats12Needle.mat','pHats12')
save('pCImu12Needle.mat','pCImu12')
save('pCIsig12Needle.mat','pCIsig12')
end

% need to bring in time to max viremia, max viremia, and viremia duration 
% can get cum