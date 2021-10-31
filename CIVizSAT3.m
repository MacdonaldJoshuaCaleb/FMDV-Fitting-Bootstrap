close all
load('SData1.mat')
% SynDataX = SData1(:,1:7);
% Syuntitled2nDataS = SData1(:,8:14);
% SynDataPV = SData1(:,15:end);
pvfitsC = readtable('pvfitsC.csv');
pvfitsC = table2array(pvfitsC);
pvfitsN = readtable('pvfitsN.csv');
pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsC.csv');
sfitsC = table2array(sfitsC);
sfitsN = readtable('sfitsN.csv');
sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsC.csv');
xfitsC = table2array(xfitsC);
xfitsN = readtable('xfitsN.csv');
xfitsN = table2array(xfitsN);
xsC = readtable('SAADataContact.csv');
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
params = readtable('FitParamsContact.csv');
params = table2array(params(:,4:10));
p = params(9,:);
   lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,ttt,x0);
   %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];


load('params9.mat');
load('params11.mat');   
load('params12.mat');
%load('params8.mat');

load('Fits9.mat');
load('Fits11.mat');
load('Fits12.mat');
%load('Fits8.mat')

temp1 = sort(Fits9);
temp1 = temp1(251:end-250,:);
temp2 = sort(Fits11);
temp2 = temp2(251:end-250,:);
temp3 = sort(Fits12);
temp3 = temp3(251:end-250,:);
%temp4 = sort(Fits8);
%temp4 = temp4(251:end-250,:);

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
plot(t,xsC(9,:),'b*',t,sesC(9,:),'r+',t,pvsC(9,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
%plot(tt,temp1(:,1:length(tt)),'b',tt,temp1(:,length(tt)+1:2*length(tt)),'r',tt,temp1(:,2*length(tt)+1:end),'g','linewidth',2)
hold off
title('ID12, SAT3')
%legend({'I(\tau) 95% CI','V(\tau) 95% CI', 'A(\tau 95% CI','log_{10}(I_{SAA})','FMDV','A_{VNT}','I(\tau) fit)','V(\tau) fit','A(\tau) fit'},'NumColumns',3,'FontSize',12)
ylim([0 14])
set(gca,'FontSize',16)
p = params(11,:);

   lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
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
plot(t,xsC(11,:),'b*',t,sesC(11,:),'r+',t,pvsC(11,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
title('ID16, SAT3')
hold off
ylim([0 14])
set(gca,'FontSize',16)

p = params(12,:);

   lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
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
plot(t,xsC(12,:),'b*',t,sesC(12,:),'r+',t,pvsC(12,:),'go','MarkerSize',20,'linewidth',2)
plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
title('ID17, SAT3')
hold off
ylim([0 14])
set(gca,'FontSize',16)

% p = params(8,:);
% 
%       lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
%     x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
%     atilde=0;
%     f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... untitled2
%                 (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
%                 (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
%    
%    [~,yy] = ode45(f,ttt,x0);
%    %yy = [yy(:,1)',yy(:,2)',yy(:,3)'];
% 
% figure
% hold on
% x2 = [tt,fliplr(tt)];
% inBetween = [temp4(1,1:length(tt)),fliplr(temp4(end,1:length(tt)))];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% inBetween2 = [temp4(1,1+length(tt):2*length(tt)),fliplr(temp4(end,1+length(tt):2*length(tt)))];
% h2 = fill(x2, inBetween2,'r','LineStyle','none');
% set(h2,'facealpha',.4)
% inBetween3 = [temp4(1,1+2*length(tt):end),fliplr(temp4(end,1+2*length(tt):end))];
% h3 = fill(x2, inBetween3,'g','LineStyle','none');
% set(h3,'facealpha',.4)
% plot(t,xsC(8,:),'b*',t,sesC(8,:),'r+',t,pvsC(8,:),'go','MarkerSize',20,'linewidth',2)
% plot(ttt,yy(:,1),'b',ttt,yy(:,2),'r',ttt,yy(:,3),'g','linewidth',2)
% title('ID29, SAT2')
% hold off
% ylim([0 14])
% set(gca,'FontSize',16)

%load('MeansSAT1C.mat')
%tempM = sort(MeansSAT1C);
%tempM = tempM(251:end-250,:);
%array2table(tempM(1,:))
%array2table(tempM(end,:))



for j = 1:10000
 RE9(j,:) = abs(params(9,:)-params9(j,:))./params(9,:);
end
ARE9 = (sum(RE9)/10000)*100;
array2table(ARE9)

for j = 1:10000
 RE11(j,:) = abs(params(11,:)-params11(j,:))./params(11,:);
end
ARE11= (sum(RE11)/10000)*100;
array2table(ARE11)

for j = 1:10000
 RE12(j,:) = abs(params(12,:)-params12(j,:))./params(12,:);
end
ARE12= (sum(RE12)/10000)*100;
array2table(ARE12)

% for j = 1:10000
%  RE8(j,:) = abs(params(8,:)-params8(j,:))./params(8,:);
% end
% ARE8= (sum(RE8)/10000)*100;
% array2table(ARE8)