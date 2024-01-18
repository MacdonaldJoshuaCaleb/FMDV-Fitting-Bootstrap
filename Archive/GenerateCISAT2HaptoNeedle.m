function [] = GenerateCISAT2HaptoNeedle
close all
%options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsN = readtable('pvfitsNHapto.csv');
pvfitsN = table2array(pvfitsN);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsN = readtable('sfitsNHapto.csv');
sfitsN = table2array(sfitsN);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsN = readtable('xfitsNHapto.csv');
xfitsN = table2array(xfitsN);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsN = readtable('HaptoDataNeedle.csv');
xsN = table2array(xsN(:,3:end));
xsN(xsN>0) = log10(xsN(xsN>0));
sesN = readtable('ViremiaNeedle.csv');
sesN = table2array(sesN(:,3:end));
pvsC = readtable('VNTNeedleData.csv');
pvsC = table2array(pvsC(:,3:end));
% xsN = readtable('SAADataNeedle.csv');
% xsN = table2array(xsN(:,3:end));
% xsN(xsN>0) = log10(xsN(xsN>0));
% sesN = readtable('ViremiaNeedle.csv');
% sesN = table2array(sesN(:,3:end));
% pvsN = readtable('VNTNeedleData.csv');
% pvsN = table2array(pvsN(:,3:end));
params = readtable('FitParamsNeedleHapto.csv');
params = table2array(params(:,4:end));
tt = 0:.1:30;
its = 10000;

Fits5 = zeros(its,3*length(tt));
Fits6 = zeros(its,3*length(tt));
Fits7 = zeros(its,3*length(tt));
Fits8 = zeros(its,3*length(tt));

params5 = zeros(its,11);
params6 = zeros(its,11);
params7 = zeros(its,11);
params8 = zeros(its,11);

%cumviralContact1 = zeros(its,1);
% cumviralContact2 = zeros(its,1);
% cumviralContact3 = zeros(its,1);
% cumviralContact4 = zeros(its,1);
% 
% maxViralContact1 = zeros(its,2);
% maxViralContact2 = zeros(its,2);
% maxViralContact3 = zeros(its,2);
% maxViralContact4 = zeros(its,2);

SData5 = zeros(its,24);
SData6 = zeros(its,24);
SData7 = zeros(its,24);
SData8 = zeros(its,24);
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
function yy = paramfun1(p,t)
    vv = ViralGrowth(dd(length(xfitsN(index,:))+1:2*length(xfitsN(index,:))));
    lambda1=(1/21.23)*p(end-2);
    
    
        lambda2 = (1/21.23)*dd(length(xfitsN(index,:)));
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
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(1:end,1)',yy(1:end,2)',yy(:,3)'];
end 

for count = 1:its
for index = 5:8
nLevel = [.4,.4,.4];
dd = GetData(xfitsN(index,:),sfitsN(index,:),pvfitsN(index,:),nLevel);

% if index <= 12    % contact
%     t = [0,2, 4, 6, 9, 12, 28];
% end
%if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
%end
if index == 5
    params5(count,:) = FittingforCIHaptoNeedle(dd,index);
    Fits5(count,:) = paramfun1(params5(count,:),tt);
    %[mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    %cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    %maxViralContact1(count,:) = [mVi tt(I)];
    SData5(count,:) = dd;
end
if index == 6
    params6(count,:) = FittingforCIHaptoNeedle(dd,index);
    Fits6(count,:) = paramfun1(params6(count,:),tt);
    %[mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    %cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    %maxViralContact2(count,:) = [mVi tt(I)];
    SData6(count,:) = dd;
end
if index == 7
    params7(count,:) = FittingforCIHaptoNeedle(dd,index);
    Fits7(count,:) = paramfun1(params7(count,:),tt);
    %[mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    %cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    %maxViralContact3(count,:) = [mVi tt(I)];
    SData7(count,:) = dd;
end
if index == 8
    params8(count,:) = FittingforCIHaptoNeedle(dd,index);
    Fits8(count,:) = paramfun1(params8(count,:),tt);
    %[mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    %cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    %maxViralContact4(count,:) = [mVi tt(I)];
    SData8(count,:) = dd;
end

end
SAT2means(count,:) = (params5(count,:)+params6(count,:)+params7(count,:)+params8(count,:))./4;
SAT2devs(count,:) = std([params5(count,:); params6(count,:); params7(count,:); params8(count,:)]);
RE5(count,:) = abs(params5(count,:)-params(5,:))./params(5,:);
RE6(count,:) = abs(params6(count,:)-params(6,:))./params(6,:);
RE7(count,:) = abs(params7(count,:)-params(7,:))./params(7,:);
RE8(count,:) = abs(params8(count,:) - params(8,:))./params(8,:); 
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


save('Fits5HaptoNeedle.mat','Fits5')
save('Fits6HaptoNeedle.mat','Fits6')
save('Fits7HaptoNeedle.mat','Fits7')
save('Fits8HaptoNeedle.mat','Fits8')

save('params5HaptoNeedle.mat','params5')
save('params6HaptoNeedle.mat','params6')
save('params7HaptoNeedle.mat','params7')
save('params8HaptoNeedle.mat','params8')

%save('cumviralContact1.mat','cumviralContact1')
%save('cumviralContact2.mat','cumviralContact2')
%save('cumviralContact3.mat','cumviralContact3')
%save('cumviralContact4.mat','cumviralContact4')

%save('maxViralContact1.mat','maxViralContact1')
%save('maxViralContact2.mat','maxViralContact2')
%save('maxViralContact3.mat','maxViralContact3')
%save('maxViralContact4.mat','maxViralContact4')
ARE5 = mean(RE5)*100;
ARE6 = mean(RE6)*100;
ARE7 = mean(RE7)*100;
ARE8 = mean(RE8)*100;
array2table(ARE5)
array2table(ARE6)
array2table(ARE7)
array2table(ARE8)
for ss = 1:11
[pHat,pCI] = lognfit(params5(:,ss));
pHats5(ss,:) = pHat;
pCImu5(ss,:) = pCI(:,1);
pCIsig5(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params6(:,ss));
pHats6(ss,:) = pHat;
pCImu6(ss,:) = pCI(:,1);
pCIsig6(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params7(:,ss));
pHats7(ss,:) = pHat;
pCImu7(ss,:) = pCI(:,1);
pCIsig7(ss,:) = pCI(:,2);
[pHat,pCI] = lognfit(params8(:,ss));
pHats8(ss,:) = pHat;
pCImu8(ss,:) = pCI(:,1);
pCIsig8(ss,:) = pCI(:,2);
end

%array2table(alphas)
%array2table(betas)
for ss = 1:11
    figure
    hold on
    %histogram(params1(:,ss))
    x = .4*min(params5(:,ss)):.0001:1.6*max(params5(:,ss));
    plot(x,lognpdf(x,pHats5(ss,1),pHats5(ss,2)),'linewidth',2)
    hold off
    xlim([min(x) max(x)])
end
% array2table(alphasSAT1)
% array2table(betasSAT1)
save('SData5HaptoNeedle.mat','SData5')
save('SData6HaptoNeedle.mat','SData6')
save('SData7HaptoNeedle.mat','SData7')
save('SData8HaptoNeedle.mat','SData8')
save('SAT2meansHaptoNeedle.mat','SAT2means')
save('SAT2devsHaptoNeedle.mat','SAT2devs')
save('pHats5HaptoNeedle.mat','pHats5')
save('pCImu5HaptoNeedle.mat','pCImu5')
save('pCIsig5HaptoNeedle.mat','pCIsig5')
save('pHats6HaptoNeedle.mat','pHats6')
save('pCImu6HaptoNeedle.mat','pCImu6')
save('pCIsig6HaptoNeedle.mat','pCIsig6')
save('pHats7HaptoNeedle.mat','pHats7')
save('pCImu7HaptoNeedle.mat','pCImu7')
save('pCIsig7HaptoNeedle.mat','pCIsig7')
save('pHats8HaptoNeedle.mat','pHats8')
save('pCImu8HaptoNeedle.mat','pCImu8')
save('pCIsig8HaptoNeedle.mat','pCIsig8')
end

% need to bring in time to max viremia, max viremia, and viremia duration 
