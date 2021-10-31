function [] = GenerateCISAT1
close all
%options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsC = readtable('pvfitsC.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsC.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsC.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('SAADataContact.csv');
xsC = table2array(xsC(:,3:end));
xsC(xsC>0) = log10(xsC(xsC>0));
sesC = readtable('ViremiaContact.csv');
sesC = table2array(sesC(:,3:end));
pvsC = readtable('VNTContactData.csv');
pvsC = table2array(pvsC(:,3:end));
xsN = readtable('SAADataNeedle.csv');
xsN = table2array(xsN(:,3:end));
xsN(xsN>0) = log10(xsN(xsN>0));
sesN = readtable('ViremiaNeedle.csv');
sesN = table2array(sesN(:,3:end));
pvsN = readtable('VNTNeedleData.csv');
pvsN = table2array(pvsN(:,3:end));
params = readtable('FitParamsContact.csv');
params = table2array(params(:,4:end));
tt = 0:.1:30;
its = 10000;

Fits1 = zeros(its,3*length(tt));
Fits2 = zeros(its,3*length(tt));
Fits3 = zeros(its,3*length(tt));
Fits4 = zeros(its,3*length(tt));

params1 = zeros(its,11);
params2 = zeros(its,11);
params3 = zeros(its,11);
params4 = zeros(its,11);

%cumviralContact1 = zeros(its,1);
% cumviralContact2 = zeros(its,1);
% cumviralContact3 = zeros(its,1);
% cumviralContact4 = zeros(its,1);
% 
% maxViralContact1 = zeros(its,2);
% maxViralContact2 = zeros(its,2);
% maxViralContact3 = zeros(its,2);
% maxViralContact4 = zeros(its,2);

SData1 = zeros(its,21);
SData2 = zeros(its,21);
SData3 = zeros(its,21);
SData4 = zeros(its,21);

function yy = paramfun1(p,t)
    vv = ViralGrowth(dd(length(xfitsC(index,:))+1:2*length(xfitsC(index,:))));
    lambda1=(1/11.18)*p(end-2);
    
    
        lambda2 = (1/11.18)*dd(length(xfitsC(index,:)));
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/11.18; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),vv(1),1e-4*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0);
   yy = [yy(1:end,1)',yy(1:end,2)',yy(:,3)'];
end 

for count = 1:its
for index = 1:4
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 1
    params1(count,:) = FittingforCI(dd,index);
    Fits1(count,:) = paramfun1(params1(count,:),tt);
    %[mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    %cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    %maxViralContact1(count,:) = [mVi tt(I)];
    SData1(count,:) = dd;
end
if index == 2
    params2(count,:) = FittingforCI(dd,index);
    Fits2(count,:) = paramfun1(params2(count,:),tt);
    %[mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    %cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    %maxViralContact2(count,:) = [mVi tt(I)];
    SData2(count,:) = dd;
end
if index == 3
    params3(count,:) = FittingforCI(dd,index);
    Fits3(count,:) = paramfun1(params3(count,:),tt);
    %[mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    %cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    %maxViralContact3(count,:) = [mVi tt(I)];
    SData3(count,:) = dd;
end
if index == 4
    params4(count,:) = FittingforCI(dd,index);
    Fits4(count,:) = paramfun1(params4(count,:),tt);
    %[mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    %cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    %maxViralContact4(count,:) = [mVi tt(I)];
    SData4(count,:) = dd;
end

end
SAT1means(count,:) = (params1(count,:)+params2(count,:)+params3(count,:)+params4(count,:))./4;
SAT1devs(count,:) = std([params1(count,:); params2(count,:); params3(count,:); params4(count,:)]);
RE1(count,:) = abs(params1(count,:)-params(1,:))./params(1,:);
RE2(count,:) = abs(params2(count,:)-params(2,:))./params(2,:);
RE3(count,:) = abs(params3(count,:)-params(3,:))./params(3,:);
RE4(count,:) = abs(params4(count,:) - params(4,:))./params(4,:); 
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
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end


save('Fits1.mat','Fits1')
save('Fits2.mat','Fits2')
save('Fits3.mat','Fits3')
save('Fits4.mat','Fits4')

save('params1.mat','params1')
save('params2.mat','params2')
save('params3.mat','params3')
save('params4.mat','params4')

%save('cumviralContact1.mat','cumviralContact1')
%save('cumviralContact2.mat','cumviralContact2')
%save('cumviralContact3.mat','cumviralContact3')
%save('cumviralContact4.mat','cumviralContact4')

%save('maxViralContact1.mat','maxViralContact1')
%save('maxViralContact2.mat','maxViralContact2')
%save('maxViralContact3.mat','maxViralContact3')
%save('maxViralContact4.mat','maxViralContact4')
ARE1 = mean(RE1)*100;
ARE2 = mean(RE2)*100;
ARE3 = mean(RE3)*100;
ARE4 = mean(RE4)*100;
array2table(ARE1)
array2table(ARE2)
array2table(ARE3)
array2table(ARE4)
for ss = 1:11
[pHat,pCI] = lognfit(params1(:,ss));
pHats1(ss,:) = pHat;
pCImu1(ss,:) = pCI(:,1);
pCIsig1(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params2(:,ss));
pHats2(ss,:) = pHat;
pCImu2(ss,:) = pCI(:,1);
pCIsig2(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params3(:,ss));
pHats3(ss,:) = pHat;
pCImu3(ss,:) = pCI(:,1);
pCIsig3(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params4(:,ss));
pHats4(ss,:) = pHat;
pCImu4(ss,:) = pCI(:,1);
pCIsig4(ss,:) = pCI(:,2);
end

%array2table(alphas)
%array2table(betas)
for ss = 1:11
    figure
    hold on
    %histogram(params1(:,ss))
    x = .4*min(params1(:,ss)):.0001:1.6*max(params1(:,ss));
    plot(x,lognpdf(x,pHats1(ss,1),pHats1(ss,2)),'linewidth',2)
    hold off
    xlim([min(x) max(x)])
end
% array2table(alphasSAT1)
% array2table(betasSAT1)
save('SData1.mat','SData1')
save('SData2.mat','SData2')
save('SData3.mat','SData3')
save('SData4.mat','SData4')
save('SAT1means.mat','SAT1means')
save('SAT1devs.mat','SAT1devs')
save('pHats1.mat','pHats1')
save('pCImu1.mat','pCImu1')
save('pCIsig1.mat','pCIsig1')
save('pHats2.mat','pHats2')
save('pCImu2.mat','pCImu2')
save('pCIsig2.mat','pCIsig2')
save('pHats3.mat','pHats3')
save('pCImu3.mat','pCImu3')
save('pCIsig3.mat','pCIsig3')
save('pHats4.mat','pHats4')
save('pCImu4.mat','pCImu4')
save('pCIsig4.mat','pCIsig4')
end

% need to bring in time to max viremia, max viremia, and viremia duration 
% can get cum