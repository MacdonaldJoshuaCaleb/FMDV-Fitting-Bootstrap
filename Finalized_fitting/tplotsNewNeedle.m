function [] = tplotsNewNeedle(its)
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


    xsC = readtable('HaptoDataNeedle.csv');
    dec = 21.23;

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

% load('paramsFinalNewAdjust.mat')
load('params_noise_needle.mat')
% load('params.mat')
params = params_noise_needle;
% params = params;
%params = paramsNoise;
    %params = readtable('FitParamsHapto.csv');

%params = params(:,4:end);
%Starts = readtable('InfectionStartTimes.csv');
%Starts = table2array(Starts);
Starts = params(:,end,its);

Start = max(Starts);
% t = [0,2, 4, 6, 9,12,28];
% t2 = t;
 t = [0,2, 4, 6, 8, 11, 14, 30];
 t2=t;
ind = find(t > Start);




x = xsC(its,:);


% 
% p = table2array(params(its,1:8));
%    lambda1=(1/dec)*p(end-3);
%     
%     
%         lambda2 = (1/dec)*x(end);
%         lambda = min(lambda1,lambda2);
%     %lambda = lambda2;
%     k=p(1); d=1/dec; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
%     x0 = [p(5),p(8),1e-4*p(5)]; % initial conditions 
%     atilde=0;
%  %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
%   %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
%    %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
%      f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
%                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
%                 (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
%    
%    [~,yy] = ode45(f,tt,x0,options);
%    base = yy;
for j = 1:10000
    tTest = [Start,t(ind)];
    tt = tTest(1):.1:20;
    tt1 = 0;
    if Starts(j) < Start
    tt1 = linspace(Starts(j),Start-(Start-Starts(j))/2,50);
    tt2 = tTest(1):.1:20;
    tt = [tt1,tt2];
    end
    p = params(j,:,its);

    
    
        lambda = (1/dec)*min(x)*.8;

    %lambda = lambda2;
    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); % use disperse here ...
    x0 = [p(7),p(9),p(8)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,tt,x0,options);
   if length(tt1) == 1
   sols1(j,:) = yy(:,1); % innate
   sols2(j,:) = yy(:,2); % pathogen
   sols3(j,:) = yy (:,3); % adaptive 
   end
   if length(tt1) > 1
       sols1(j,:) = yy(length(tt1)+1:end,1); % innate
       sols2(j,:) = yy(length(tt1)+1:end,2); % pathogen
       sols3(j,:) = yy (length(tt1)+1:end,3); % adaptive 
   end
 
  
   if mod(j,100) == 0
       fprintf('iteration %i of %i\n',j,10000)
   end
end
remove = 250;

sols1 = sort(sols1);
sols2 = sort(sols2);
sols3 = sort(sols3);
 %gam = sort(gam);
% bet = sort(bet);
lb1 = sols1(remove+1,:);
ub1 = sols1(end-remove,:);
lb2 = sols2(remove+1,:);
ub2 = sols2(end-remove,:);
lb3 = sols3(remove+1,:);
ub3 = sols3(end-remove,:);
% lb4 = gam(remove+1,:);
% ub4 = gam(end-remove,:);
% lb5 = bet(remove+1,:);
% ub5 = bet(end-remove,:);
base(:,1) = median(sols1);
base(:,2) = median(sols2);
base(:,3) = median(sols3);
% base(:,4) = mean(gam);
% base(:,5) = mean(bet);
str = strcat('ID',' ',num2str(IDs(its)),{', '},extractBetween(groups(its),1,4));
    str2= strcat('ID',num2str(IDs(its)),groups(its),'CIHaptoNeedle');

   %fig = figure;
   red = [0.9290, 0.6940, 0.1250];
   blue = [0 0 0];
   left_color = blue;
   right_color = red;
   %set(fig,'defaultAxesColorOrder',[left_color; right_color]);
%plot(tt,lb1,'b',tt,ub1,'b','linewidth',2)
%yyaxis left
tt = Start:.1:20;
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

%     fever = readtable('QualCompInd.csv');
%     fstart = table2array(fever(its,end-2))-2;
%     fend = table2array(fever(its,end-1))-2;
%     if fstart > 0
%     xf = [fstart,fend];
%     xf = [xf,fliplr(xf)];
%     inBetween = [[0,0],fliplr([11,11])];
%     h = fill(xf, inBetween,'y','LineStyle','none','marker','none');
%     set(h,'facealpha',.2)
%     %text(fstart,9,"fever detected")
%     end
% xline(Start,'k--')
% y = randsample(1:10000,20);
% for j = 1:length(y)
%     tTest = [Start,t(ind)];
%     tt = tTest(1):.1:15;
%     tt1 = 0;
%     if Starts(y(j)) < Start
%     tt1 = linspace(Starts(y(j)),Start,50);
%        p = params(y(j),:,its);
% 
%     
%     
%         lambda = (1/dec)*min(x).*.8;
% 
%     %lambda = lambda2;
%     k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); % use disperse here ...
%     x0 = [p(7),p(9),p(8)]; % initial conditions 
%     atilde=0;
%  %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
%   %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
%    %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
%      f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
%                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
%                 (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
%    
%    [~,yy] = ode45(f,tt1,x0,options);
%     end
%     plot(tt1,yy(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
%     plot(tt1,yy(:,2),'Color',cmap(3,:),'LineStyle',':','linewidth',2)
%     plot(tt1,yy(:,3),'Color',cmap(2,:),'LineStyle',':','linewidth',2)
% end
hold off
ylabel('Concentration')
ylim([0 11])
% yyaxis right
% plot(TempTimes,Temps(:,its),'color',[0.9290, 0.6940, 0.1250],'linewidth',2)
 xlim([0 20])
% ylim([-2.05, 2.75])
% ylabel('Temprature residuals')
xlabel('Time since contact was initiated (days)')
title(str)
set(gca,'FontSize',16)
baseFileName = sprintf(str2);
fname = '~/Documents/MATLAB/Model_identifiability/Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');


end