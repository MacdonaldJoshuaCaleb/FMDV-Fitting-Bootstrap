function [] = GenerateCISAT3HAPTO
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

Fits9 = zeros(its,3*length(tt));
%Fits10 = zeros(its,3*length(tt));
Fits11 = zeros(its,3*length(tt));
Fits12 = zeros(its,3*length(tt));

params9 = zeros(its,7);
% params10 = zeros(its,7);
params11 = zeros(its,7);
params12 = zeros(its,7);

cumviralContact9 = zeros(its,1);
% cumviralContact10 = zeros(its,1);
cumviralContact11 = zeros(its,1);
cumviralContact12 = zeros(its,1);

maxViralContact9 = zeros(its,2);
% maxViralContact10 = zeros(its,2);
maxViralContact11 = zeros(its,2);
maxViralContact12 = zeros(its,2);

SData9 = zeros(its,21);
% SData10 = zeros(its,21);
SData11 = zeros(its,21);
SData12 = zeros(its,21);

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
for index = 9:12
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 9
    params9(count,:) = FittingforCIHapto(dd,index);
    Fits9(count,:) = paramfun1(params9(count,:),tt);
    [mVi I] = max(Fits9(length(tt)+1:2*length(tt)));
    cumviralContact9(count) = sum(Fits9(count,length(tt)+1:2*length(tt)));
    maxViralContact9(count,:) = [mVi tt(I)];
    SData9(count,:) = dd;
end
% if index == 10
%     params10(count,:) = FittingforCI(dd,index);
%     Fits10(count,:) = paramfun1(params10(count,:),tt);
%     [mVi I] = max(Fits10(length(tt)+1:2*length(tt)));
%     cumviralContact10(count) = sum(Fits10(count,length(tt)+1:2*length(tt)));
%     maxViralContact10(count,:) = [mVi tt(I)];
%     SData10(count,:) = dd;
% end
if index == 11
    params11(count,:) = FittingforCIHapto(dd,index);
    Fits11(count,:) = paramfun1(params11(count,:),tt);
    [mVi I] = max(Fits11(length(tt)+1:2*length(tt)));
    cumviralContact11(count) = sum(Fits11(count,length(tt)+1:2*length(tt)));
    maxViralContact11(count,:) = [mVi tt(I)];
    SData11(count,:) = dd;
end
if index == 12
    params12(count,:) = FittingforCIHapto(dd,index);
    Fits12(count,:) = paramfun1(params12(count,:),tt);
    [mVi I] = max(Fits12(length(tt)+1:2*length(tt)));
    cumviralContact12(count) = sum(Fits12(count,length(tt)+1:2*length(tt)));
    maxViralContact12(count,:) = [mVi tt(I)];
    SData12(count,:) = dd;
end

end
fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end



save('Fits9Hapto.mat','Fits9')
%save('Fits10.mat','Fits10')
save('Fits11Hapto.mat','Fits11')
save('Fits12Hapto.mat','Fits12')

save('params9Hapto.mat','params9')
%save('params10.mat','params10')
save('params11Hapto.mat','params11')
save('params12Hapto.mat','params12')

save('cumviralContact9Hapto.mat','cumviralContact9')
%save('cumviralContact10.mat','cumviralContact10')
save('cumviralContact11Hapto.mat','cumviralContact11')
save('cumviralContact12Hapto.mat','cumviralContact12')

save('maxViralContact9Hapto.mat','maxViralContact9')
%save('maxViralContact10.mat','maxViralContact10')
save('maxViralContact11Hapto.mat','maxViralContact11')
save('maxViralContact12Hapto.mat','maxViralContact12')

save('SData9Hapto.mat','SData9')
%save('SData10.mat','SData10')
save('SData11Hapto.mat','SData11')
save('SData12.matHapto','SData12')
end