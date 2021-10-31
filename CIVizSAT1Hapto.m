close all
load('SData1.mat')
% SynDataX = SData1(:,1:7);
% Syuntitled2nDataS = SData1(:,8:14);
% SynDataPV = SData1(:,15:end);
pvfitsC = readtable('pvfitsCHapto.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsNHapto.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsCHapto.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsNHapto.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsCHapto.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('HaptoDataContact.csv');
xsC = table2array(xsC(:,3:end));
xsC(xsC>0) = log10(xsC(xsC>0));
sesC = readtable('ViremiaContact.csv');   load('SData1.mat')
% SynDataX = SData1(:,1:7);
% SynDataS = SData1(:,8:14);
% SynDataPV = SData1(:,15:end);
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
t = [0,2, 4, 6, 9, 12, 28];
tt = 0:.1:30;
ttt=0:.01:30;
params = readtable('FitParamsHapto.csv');
params = table2array(params(:,:));
p = params(1,:);
   lambda=0; k=p(3); d=1/21.23; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,ttt,x0);
   %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];


load('params1Hapto.mat');
load('params2Hapto.mat');   load('SData1.mat')
SynDataX = SData1(:,1:7);
SynDataS = SData1(:,8:14);
SynDataPV = SData1(:,15:end);
load('params4Hapto.mat');
load('params3Hapto.mat');

load('Fits1Hapto.mat');
load('Fits2Hapto.mat');
load('Fits3Hapto.mat');
load('Fits4Hapto.mat')

temp1 = sort(Fits1);
temp1 = temp1(251:end-250,:);
temp2 = sort(Fits2);
temp2 = temp2(251:end-250,:);
temp3 = sort(Fits3);
temp3 = temp3(251:end-250,:);
temp4 = sort(Fits4);
temp4 = temp4(251:end-250,:);

figure
hold on
x2 = [tt,fliplr(tt)];
inBetween = [temp1(1,1:length(tt)),fliplr(temp1(end,1:length(tt)))];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
inBetween2 = [temp1(1,1+length(tt):2*length(tt)),fliplr(temp1(end,1+length(tt):2*length(tt)))];
h2 = fill(x2, inBetween2,'r','LineStyle','none');
set(h2,'facealpha',.4)
inBetween3 = [temp1(1,1+2*length(tt):end),fliplr(temp1(end,1+2*length(tt):end))];
h3 = fill(x2, inBetween3,'g','LineStyle','none');
set(h3,'facealpha',.4)
plot(t,xsC(1,:),'b*',t,sesC(1,:),'r+',t,pvsC(1,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
%plot(tt,temp1(:,1:length(tt)),'b',tt,temp1(:,length(tt)+1:2*length(tt)),'r',tt,temp1(:,2*length(tt)+1:end),'g','linewidth',2)
hold off
title('ID2, SAT1')
legend({'I(\tau) 95% CI','V(\tau) 95% CI', 'A(\tau 95% CI','log_{10}(I_{Hapto})','FMDV','A_{VNT}','I(\tau) fit)','V(\tau) fit','A(\tau) fit'},'NumColumns',3,'FontSize',12)
ylim([0 14])
set(gca,'FontSize',16)
p = params(2,:);

   lambda=0; k=p(3); d=1/21.23; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,ttt,x0);
   %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];

figure
hold on
x2 = [tt,fliplr(tt)];
inBetween = [temp2(1,1:length(tt)),fliplr(temp2(end,1:length(tt)))];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
inBetween2 = [temp2(1,1+length(tt):2*length(tt)),fliplr(temp2(end,1+length(tt):2*length(tt)))];
h2 = fill(x2, inBetween2,'r','LineStyle','none');
set(h2,'facealpha',.4)
inBetween3 = [temp2(1,1+2*length(tt):end),fliplr(temp2(end,1+2*length(tt):end))];
h3 = fill(x2, inBetween3,'g','LineStyle','none');
set(h3,'facealpha',.4)
plot(t,xsC(2,:),'b*',t,sesC(2,:),'r+',t,pvsC(2,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
title('ID4, SAT1')
hold off
ylim([0 14])
set(gca,'FontSize',16)

p = params(3,:);

     lambda=0; k=p(3); d=1/21.23; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   [~,yy] = ode45(f,ttt,x0);
   %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];

figure
hold on
x2 = [tt,fliplr(tt)];
inBetween = [temp3(1,1:length(tt)),fliplr(temp3(end,1:length(tt)))];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
inBetween2 = [temp3(1,1+length(tt):2*length(tt)),fliplr(temp3(end,1+length(tt):2*length(tt)))];
h2 = fill(x2, inBetween2,'r','LineStyle','none');
set(h2,'facealpha',.4)
inBetween3 = [temp3(1,1+2*length(tt):end),fliplr(temp3(end,1+2*length(tt):end))];
h3 = fill(x2, inBetween3,'g','LineStyle','none');
set(h3,'facealpha',.4)
plot(t,xsC(3,:),'b*',t,sesC(3,:),'r+',t,pvsC(3,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
title('ID19, SAT1')
hold off
ylim([0 14])
set(gca,'FontSize',16)

p = params(4,:);

      lambda=0; k=p(3); d=1/21.23; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... untitled2
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,ttt,x0);
   %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];

figure
hold on
x2 = [tt,fliplr(tt)];
inBetween = [temp4(1,1:length(tt)),fliplr(temp4(end,1:length(tt)))];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
inBetween2 = [temp4(1,1+length(tt):2*length(tt)),fliplr(temp4(end,1+length(tt):2*length(tt)))];
h2 = fill(x2, inBetween2,'r','LineStyle','none');
set(h2,'facealpha',.4)
inBetween3 = [temp4(1,1+2*length(tt):end),fliplr(temp4(end,1+2*length(tt):end))];
h3 = fill(x2, inBetween3,'g','LineStyle','none');
set(h3,'facealpha',.4)
plot(t,xsC(4,:),'b*',t,sesC(4,:),'r+',t,pvsC(4,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
title('ID33, SAT1')
hold off
ylim([0 14])
set(gca,'FontSize',16)

load('MeansSAT1C.mat')
tempM = sort(MeansSAT1C);
tempM = tempM(251:end-250,:);
array2table(tempM(1,:))
array2table(tempM(end,:))




for j = 1:10000
 RE1(j,:) = abs(params(1,1:end-1)-params1(j,:))./params(1,1:end-1);
end
ARE1 = (sum(RE1)/10000)*100;
array2table(ARE1)

for j = 1:10000
 RE2(j,:) = abs(params(2,1:end-1)-params2(j,:))./params(2,1:end-1);
end
ARE2= (sum(RE2)/10000)*100;
array2table(ARE2)

for j = 1:10000
 RE3(j,:) = abs(params(3,1:end-1)-params3(j,:))./params(3,1:end-1);
end
ARE3= (sum(RE3)/10000)*100;
array2table(ARE3)

for j = 1:10000
 RE4(j,:) = abs(params(4,1:end-1)-params4(j,:))./params(4,1:end-1);
end
ARE4= (sum(RE4)/10000)*100;
array2table(ARE4)