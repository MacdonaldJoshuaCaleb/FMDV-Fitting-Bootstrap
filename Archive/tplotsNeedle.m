function [] = tplotsNeedle(its,opt)
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
xsC = readtable('SAADataNeedle.csv');
dec = 11.18;
end
if opt == 2
    xsC = readtable('HaptoDataNeedle.csv');
    dec = 21.23;
end
groups = xsC(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xsC(:,2);
IDs = table2array(IDs);
xsC = table2array(xsC(:,3:end));
xsC(xsC>0) = log10(xsC(xsC>0));
sesC = readtable('ViremiaNeedle.csv');  
sesC = table2array(sesC(:,3:end));
pvsC = readtable('VNTNeedleData.csv');
pvsC = table2array(pvsC(:,3:end));
if opt == 1
load('params1Needle.mat')
load('params2Needle.mat')
load('params3Needle.mat')
load('params4Needle.mat')
load('params5Needle.mat')
load('params6Needle.mat')
load('params7Needle.mat')
load('params8Needle.mat')
load('params9Needle.mat')
load('params10Needle.mat')
load('params11Needle.mat')
load('params12Needle.mat')
params = readtable('FitParamsNeedle.csv');
end
if opt == 2
    load('params1HaptoNeedle.mat')
load('params2HaptoNeedle.mat')
load('params3HaptoNeedle.mat')
load('params4HaptoNeedle.mat')
load('params5HaptoNeedle.mat')
load('params6HaptoNeedle.mat')
load('params7HaptoNeedle.mat')
load('params8HaptoNeedle.mat')
load('params9HaptoNeedle.mat')
load('params10HaptoNeedle.mat')
load('params11HaptoNeedle.mat')
load('params12HaptoNeedle.mat')
    params = readtable('FitParamsNeedleHapto.csv');
end
params = params(:,4:end);
%t = [0,2, 4, 6, 9,12,28];
 t = [0,2, 4, 6, 8, 11, 14, 30];
t2 = t;
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
if its == 10
    ps = params10(:,1:8);
end
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
   sols1(j,:) = yy(:,1);
   sols2(j,:) = yy(:,2);
   sols3(j,:) = yy (:,3);
   if mod(j,100) == 0
       fprintf('iteration %i of %i\n',j,10000)
   end
end
remove = 250;
sols1 = sort(sols1);
sols2 = sort(sols2);
sols3 = sort(sols3);
lb1 = sols1(remove+1,:);
ub1 = sols1(end-remove,:);
lb2 = sols2(remove+1,:);
ub2 = sols2(end-remove,:);
lb3 = sols3(remove+1,:);
ub3 = sols3(end-remove,:);
base(:,1) = mean(sols1);
base(:,2) = mean(sols2);
base(:,3) = mean(sols3);
str = strcat('ID',' ',num2str(IDs(its)),', ',extractBetween(groups(its),1,4));
if opt == 1
str2= strcat('ID',num2str(IDs(its)),groups(its),'CISAA');
end
if opt == 2
    str2= strcat('ID',num2str(IDs(its)),groups(its),'CIHapto');
end
figure 
hold on
%plot(tt,lb1,'b',tt,ub1,'b','linewidth',2)
hold on
x2 = [tt,fliplr(tt)];
inBetween = [lb1,fliplr(ub1)];
h = fill(x2, inBetween,cmap(1,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
scatter(t,x,120,cmap(1,:),'Marker','+','linewidth',2)
%plot(tt,lb2,'r',tt,ub2,'r','linewidth',2)
inBetween = [lb2,fliplr(ub2)];
h = fill(x2, inBetween,cmap(3,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,2),'Color',cmap(3,:),'LineStyle','-','linewidth',2)
scatter(t,sesC(its,:),120,cmap(3,:),'Marker','*','linewidth',2)
%plot(tt,lb3,'g',tt,ub3,'g','linewidth',2)
inBetween = [lb3,fliplr(ub3)];
h = fill(x2, inBetween,cmap(2,:),'LineStyle','none');
set(h,'facealpha',.4)
scatter(t,pvsC(its,:),120,cmap(2,:),'Marker','o','linewidth',2)
plot(tt,base(:,3),'Color',cmap(2,:),'LineStyle','--','linewidth',2)
hold off
ylim([0 11])
xlim([0 15])
ylabel('Concentration')
xlabel('Days')
title(str)
set(gca,'FontSize',16)
baseFileName = sprintf(str2);
fname = '~/Documents/MATLAB/Model_identifiability/Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
end