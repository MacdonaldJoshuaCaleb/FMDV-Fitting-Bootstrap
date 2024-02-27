function [] = tplotsNewNeedle(its)
% input specify host 1-12
% output time plots in figure S3

options=odeset('RelTol',1e-12,'AbsTol',1e-12);
% define colormaps
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

% read in needle haptoglibin data 
    xsC = readtable('HaptoDataNeedle.csv');
    dec = 21.23;

% grab out serotype and id info for plot titles 
groups = xsC(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xsC(:,2);
IDs = table2array(IDs);
xsC = table2array(xsC(:,3:end));
% take log 10 of positive values
xsC(xsC>0) = log10(xsC(xsC>0));

% read in needle virus data 
sesC = readtable('ViremiaNeedle.csv');  
sesC = table2array(sesC(:,3:end));

% read in needle VNT data 
pvsC = readtable('VNTNeedleData.csv');
pvsC = table2array(pvsC(:,3:end));

% read in the parameters 
load('params_noise_needle.mat')


params = params_noise_needle;

% find start times 
Starts = params(:,end,its);

Start = max(Starts);

 t = [0,2, 4, 6, 8, 11, 14, 30];
 t2=t;
ind = find(t > Start);




x = xsC(its,:);



for j = 1:10000
    % set up time stepping so times align and CIs can be plotted 
    tTest = [Start,t(ind)];
    tt = tTest(1):.1:20;
    tt1 = 0;
    if Starts(j) < Start
    tt1 = linspace(Starts(j),Start-(Start-Starts(j))/2,50);
    tt2 = tTest(1):.1:20;
    tt = [tt1,tt2];
    end
    % solve the ode
    p = params(j,:,its);

    
    
        lambda = (1/dec)*min(x)*.8;


    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); 
    x0 = [p(7),p(9),p(8)]; % initial conditions 
    atilde=0;

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
 
  % print progress 
   if mod(j,100) == 0
       fprintf('iteration %i of %i\n',j,10000)
   end
end
remove = 250;

% Get CIs
sols1 = sort(sols1);
sols2 = sort(sols2);
sols3 = sort(sols3);

% CI LBs and UBs for each compartment 
lb1 = sols1(remove+1,:);
ub1 = sols1(end-remove,:);
lb2 = sols2(remove+1,:);
ub2 = sols2(end-remove,:);
lb3 = sols3(remove+1,:);
ub3 = sols3(end-remove,:);

base(:,1) = median(sols1);
base(:,2) = median(sols2);
base(:,3) = median(sols3);

% string names for file save
str = strcat('ID',' ',num2str(IDs(its)),{', '},extractBetween(groups(its),1,4));
    str2= strcat('ID',num2str(IDs(its)),groups(its),'CIHaptoNeedle');

% time stepping after max infection start time 
tt = Start:.1:20;

figure
hold on
% plot innate 
x2 = [tt,fliplr(tt)];
inBetween = [lb1,fliplr(ub1)];
h = fill(x2, inBetween,cmap(1,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
scatter(t2,x,120,cmap(1,:),'Marker','+','linewidth',2)

% plot virus 
inBetween = [lb2,fliplr(ub2)];
h = fill(x2, inBetween,cmap(3,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,2),'Color',cmap(3,:),'LineStyle','-','linewidth',2)
scatter(t2,sesC(its,:),120,cmap(3,:),'Marker','*','linewidth',2)

% plot adaptive 
inBetween = [lb3,fliplr(ub3)];
h = fill(x2, inBetween,cmap(2,:),'LineStyle','none','marker','none');
set(h,'facealpha',.4)
scatter(t2,pvsC(its,:),120,cmap(2,:),'Marker','o','linewidth',2)
plot(tt,base(:,3),'Color',cmap(2,:),'LineStyle','--','marker','none','linewidth',2)

% plot formatting 
hold off
ylabel('Concentration')
ylim([0 11])
xlim([0 20])

xlabel('Time since contact was initiated (days)')
title(str)
set(gca,'FontSize',16)

% save the figure 
baseFileName = sprintf(str2);
fname = '~/Documents/MATLAB/Model_identifiability/Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');


end