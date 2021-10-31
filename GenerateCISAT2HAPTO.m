function [] = GenerateCISAT2HAPTO
close all
%options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsC = readtable('pvfitsCHapto.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsCHapto.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsCHapto.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('HaptoDataContact.csv');
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

tt = 0:.1:30;
its = 10000;

Fits5 = zeros(its,3*length(tt));
Fits6 = zeros(its,3*length(tt));
Fits7 = zeros(its,3*length(tt));
Fits8 = zeros(its,3*length(tt));

params5 = zeros(its,7);
params6 = zeros(its,7);
params7 = zeros(its,7);
params8 = zeros(its,7);

cumviralContact5 = zeros(its,1);
cumviralContact6 = zeros(its,1);
cumviralContact7 = zeros(its,1);
cumviralContact8 = zeros(its,1);

maxViralContact5 = zeros(its,2);
maxViralContact6 = zeros(its,2);
maxViralContact7 = zeros(its,2);
maxViralContact8 = zeros(its,2);

SData5 = zeros(its,21);
SData6 = zeros(its,21);
SData7 = zeros(its,21);
SData8 = zeros(its,21);

function yy = paramfun1(p,t)
         lambda=0; k=p(3); d=1/21.23; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0);
   yy = [yy(:,1)',yy(:,2)',yy(:,3)'];
end 

for count = 1:its
for index = 5:8
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 5
    params5(count,:) = FittingforCIHapto(dd,index);
    Fits5(count,:) = paramfun1(params5(count,:),tt);
    [mVi I] = max(Fits5(length(tt)+1:2*length(tt)));
    cumviralContact5(count) = sum(Fits5(count,length(tt)+1:2*length(tt)));
    maxViralContact5(count,:) = [mVi tt(I)];
    SData5(count,:) = dd;
end
if index == 6
    params6(count,:) = FittingforCIHapto(dd,index);
    Fits6(count,:) = paramfun1(params6(count,:),tt);
    [mVi I] = max(Fits6(length(tt)+1:2*length(tt)));
    cumviralContact6(count) = sum(Fits6(count,length(tt)+1:2*length(tt)));
    maxViralContact6(count,:) = [mVi tt(I)];
    SData6(count,:) = dd;
end
if index == 7
    params7(count,:) = FittingforCIHapto(dd,index);
    Fits7(count,:) = paramfun1(params7(count,:),tt);
    [mVi I] = max(Fits7(length(tt)+1:2*length(tt)));
    cumviralContact7(count) = sum(Fits7(count,length(tt)+1:2*length(tt)));
    maxViralContact7(count,:) = [mVi tt(I)];
    SData7(count,:) = dd;
end
% if index == 8
%     params8(count,:) = FittingforCIHapto(dd,index);
%     Fits8(count,:) = paramfun1(params8(count,:),tt);
%     [mVi I] = max(Fits8(length(tt)+1:2*length(tt)));
%     cumviralContact8(count) = sum(Fits8(count,length(tt)+1:2*length(tt)));
%     maxViralContact8(count,:) = [mVi tt(I)];
%     SData8(count,:) = dd;
% end

end
fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end



save('Fits5Hapto.mat','Fits5')
save('Fits6Hapto.mat','Fits6')
save('Fits7Hapto.mat','Fits7')
% save('Fits8Hapto.mat','Fits8')

save('params5Hapto.mat','params5')
save('params6Hapto.mat','params6')
save('params7Hapto.mat','params7')
% save('params8Hapto.mat','params8')

save('cumviralContact5Hapto.mat','cumviralContact5')
save('cumviralContact6Hapto.mat','cumviralContact6')
save('cumviralContact7Hapto.mat','cumviralContact7')
% save('cumviralContact8.mat','cumviralContact8')

save('maxViralContact5Hapto.mat','maxViralContact5')
save('maxViralContact6Hapto.mat','maxViralContact6')
save('maxViralContact7Hapto.mat','maxViralContact7')
% save('maxViralContact8Hapto.mat','maxViralContact8')

save('SData5Hapto.mat','SData5')
save('SData6Hapto.mat','SData6')
save('SData7Hapto.mat','SData7')
% save('SData8Hapto.mat','SData8')
end