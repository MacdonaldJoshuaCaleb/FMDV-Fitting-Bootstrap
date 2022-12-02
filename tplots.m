function [] = tplots(its,opt)
%close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
% pvsfit = table2array(pvsfit);
% sesfit = table2array(sesfit);
% xsfit = table2array(xsfit);
cmap1 = [.5 .5 .5;
0 0 0;
0 0 1];
cmap2 = [.5 .5 .5;
0 0 0;
0 1 0];
cmap3 = [.5 .5 .5;
0 0 0;
1 0 0];
if its <= 4
    cmap = cmap1;
end
if its >4 && its <=8
    cmap = cmap2;
end
if its >8
    cmap = cmap3;
end

if opt == 1
xsC = readtable('SAADataContact.csv');
dec = 11.18;
end
if opt == 2
    xsC = readtable('HaptoDataContact.csv');
    dec = 21.23;
end
groups = xsC(:,1);
groups = table2array(groups);
groups = string(groups);
TData = readtable('MeanTemps.csv');
TempTimes = table2array(TData(:,1));
Temps = table2array(TData(:,2:end));
IDs = xsC(:,2);
IDs = table2array(IDs);
xsC = table2array(xsC(:,3:end));
xsC(xsC>0) = log10(xsC(xsC>0));
sesC = readtable('ViremiaContact.csv');  
sesC = table2array(sesC(:,3:end));
pvsC = readtable('VNTContactData.csv');
pvsC = table2array(pvsC(:,3:end));
if opt == 1
load('params1.mat')
load('params2.mat')
load('params3.mat')
load('params4.mat')
load('params5.mat')
load('params6.mat')
load('params7.mat')
load('params8.mat')
load('params9.mat')
%load('params10.mat')
load('params11.mat')
load('params12.mat')
params = readtable('FitParamsContact.csv');
end
if opt == 2
    load('params1Hapto.mat')
load('params2Hapto.mat')
load('params3Hapto.mat')
load('params4Hapto.mat')
load('params5Hapto.mat')
load('params6Hapto.mat')
load('params7Hapto.mat')
load('params8Hapto.mat')
load('params9Hapto.mat')
%load('params10.mat')
load('params11Hapto.mat')
load('params12Hapto.mat')
    params = readtable('FitParamsHapto.csv');
end
params = params(:,4:end);
Starts = readtable('InfectionStartTimes.csv');
Starts = table2array(Starts);
Starts = mean(Starts);
t2 = [0,2, 4, 6, 9,12,28];
 %t = [0,2, 4, 6, 8, 11, 14, 30];
t = [Starts(its),t2(2:end)];
tt = t(1):.1:30;

x = xsC(its,:);

if its == 1
    ps = params1(:,1:8);
end
if its == 2
    ps = params2(:,1:8);
end   
if  its == 3
    ps = params3(:,1:8);
end
if its == 4
    ps = params4(:,1:8);
end   
if its == 5
    ps = params5(:,1:8);
end
if its == 6
    ps = params6(:,1:8);
end
if its == 7
    ps = params7(:,1:8);
end
if its == 8
    ps = params8(:,1:8);
end
if its == 9
    ps = params9(:,1:8);
end
% if its == 10
%     ps = params10(:,1:8);
% end
if its == 11
    ps = params11(:,1:8);
end
if its == 12
    ps = params12(:,1:8);
end

p = table2array(params(its,1:8));
   lambda1=(1/dec)*p(end-3);
    
    
        lambda2 = (1/dec)*x(end);
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/dec; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),p(8),1e-4*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,tt,x0,options);
   base = yy;
for j = 1:10000
    p = ps(j,:);
    lambda1=(1/dec)*p(end-3);
    
    
        lambda2 = (1/dec)*x(end);
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/dec; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),p(8),1e-4*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,tt,x0,options);
   sols1(j,:) = yy(:,1); % innate
   sols2(j,:) = yy(:,2); % pathogen
   sols3(j,:) = yy (:,3); % adaptive 
 
   gam(j,:) = sols2(j,:).^2./(ones(1,length(sols2(j,:)))+sols2(j,:).^2);
   bet(j,:) = sols3(j,:)./(ones(1,length(sols2(j,:)))+sols2(j,:));
   if mod(j,100) == 0
       fprintf('iteration %i of %i\n',j,10000)
   end
end
remove = 250;
sols1 = sort(sols1);
sols2 = sort(sols2);
sols3 = sort(sols3);
gam = sort(gam);
bet = sort(bet);
lb1 = sols1(remove+1,:);
ub1 = sols1(end-remove,:);
lb2 = sols2(remove+1,:);
ub2 = sols2(end-remove,:);
lb3 = sols3(remove+1,:);
ub3 = sols3(end-remove,:);
lb4 = gam(remove+1,:);
ub4 = gam(end-remove,:);
lb5 = bet(remove+1,:);
ub5 = bet(end-remove,:);
base(:,1) = mean(sols1);
base(:,2) = mean(sols2);
base(:,3) = mean(sols3);
base(:,4) = mean(gam);
base(:,5) = mean(bet);
str = strcat('ID',' ',num2str(IDs(its)),', ',extractBetween(groups(its),1,4));
if opt == 1
str2= strcat('ID',num2str(IDs(its)),groups(its),'CISAA');
end
if opt == 2
    str2= strcat('ID',num2str(IDs(its)),groups(its),'CIHapto');
end
   %fig = figure;
   red = [0.9290, 0.6940, 0.1250];
   blue = [0 0 0];
   left_color = blue;
   right_color = red;
   %set(fig,'defaultAxesColorOrder',[left_color; right_color]);
%plot(tt,lb1,'b',tt,ub1,'b','linewidth',2)
%yyaxis left
figure
hold on
x2 = [tt,fliplr(tt)];
inBetween = [lb1,fliplr(ub1)];
h = fill(x2, inBetween,cmap(1,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
scatter(t2,x,120,cmap(1,:),'Marker','+','linewidth',2)
%plot(tt,lb2,'r',tt,ub2,'r','linewidth',2)
inBetween = [lb2,fliplr(ub2)];
h = fill(x2, inBetween,cmap(3,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,2),'Color',cmap(3,:),'LineStyle','-','linewidth',2)
scatter(t2,sesC(its,:),120,cmap(3,:),'Marker','*','linewidth',2)
%plot(tt,lb3,'g',tt,ub3,'g','linewidth',2)
inBetween = [lb3,fliplr(ub3)];
h = fill(x2, inBetween,cmap(2,:),'LineStyle','none','marker','none');
set(h,'facealpha',.4)
scatter(t2,pvsC(its,:),120,cmap(2,:),'Marker','o','linewidth',2)
plot(tt,base(:,3),'Color',cmap(2,:),'LineStyle','--','marker','none','linewidth',2)

    fever = readtable('QualCompInd.csv');
    fstart = table2array(fever(its,15))-2;
    fend = table2array(fever(its,16))-2;
    if fstart > 0
    xf = [fstart,fend];
    xf = [xf,fliplr(xf)];
    inBetween = [[0,0],fliplr([11,11])];
    h = fill(xf, inBetween,'y','LineStyle','none','marker','none');
    set(h,'facealpha',.2)
    %text(fstart,9,"fever detected")
    end
xline(Starts(its),'k--','linewidth',2)
hold off
ylabel('Concentration')
ylim([0 11])
% yyaxis right
% plot(TempTimes,Temps(:,its),'color',[0.9290, 0.6940, 0.1250],'linewidth',2)
 xlim([0 15])
% ylim([-2.05, 2.75])
% ylabel('Temprature residuals')
xlabel('Time since contact (days)')
title(str)
set(gca,'FontSize',16)
baseFileName = sprintf(str2);
fname = '~/Documents/MATLAB/Model_identifiability/Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');

% figure 
% hold on
% x2 = [tt,fliplr(tt)];
% inBetween = [lb4,fliplr(ub4)];
% h = fill(x2, inBetween,'r','LineStyle','none');
% set(h,'facealpha',.4)
% plot(tt,base(:,4),'r','LineStyle',':','linewidth',2)
% 
% x2 = [tt,fliplr(tt)];
% inBetween = [lb5,fliplr(ub5)];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% plot(tt,base(:,5),'b','LineStyle','--','linewidth',2)
% hold off
% 
% Q = trapz(tt,(1./tt(end)).*base(:,4));
% pfit = 3.57/Q;
% mv = [3.57, 2.33, 1.92];
% 
%    
% figure 
% hold on
% x2 = [tt,fliplr(tt)];
% inBetween = [pfit.*lb4,fliplr(pfit.*ub4)];
% h = fill(x2, inBetween,'r','LineStyle','none');
% set(h,'facealpha',.4)
% plot(tt,pfit.*base(:,4),'r','LineStyle',':','linewidth',2)
end