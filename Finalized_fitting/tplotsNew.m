function [] = tplotsNew(its,opt)
% function for plotting results for contact infected hosts 
% inputs 
% its  - specifys host 1-12, see data table
% opt - specifys plotting option 
% option 1- plots with CIs as they appear in the manuscript 
% option 2 - plots with SRS of 500 time trajectories 

options=odeset('RelTol',1e-12,'AbsTol',1e-12);
% specify serotype specific colormpas 
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
    % sat 1
    cmap = cmap1;
end
if its >4 && its <=8
    % sat 2
    cmap = cmap2;
end
if its >8
    % sat 3
    cmap = cmap3;
end

    % read in haptoglobin data
    xsC = readtable('HaptoDataContact.csv');
    % specify haptoglobin decay rate from Glidden et al.
    dec = 21.23;

    % extract host ID and serotype information 
groups = xsC(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xsC(:,2);
IDs = table2array(IDs);
xsC = table2array(xsC(:,3:end));
% log_10 transform the positive haptoglibin data values 
xsC(xsC>0) = log10(xsC(xsC>0));

% read in the qPCR virus data
sesC = readtable('ViremiaContact.csv');  
sesC = table2array(sesC(:,3:end));
% read in the VNT data 
pvsC = readtable('VNTContactData.csv');
pvsC = table2array(pvsC(:,3:end));

% read in the paramter estimates from Monte Carlo simulations 
load('params_noise.mat')

% give them a shorter name
params = params_noise;

% pull out the infection start times
Starts = params(:,end,its);

% identify the maximum start time
Start = max(Starts);
% observation time points 
t = [0,2, 4, 6, 9,12,28];
t2 = t;
 
ind = find(t > Start);



% specify the row of the data to plot 
x = xsC(its,:);
s = sesC(its,:);
pv = pvsC(its,:);

% plot display option 1 (used in manuscript) plots with CIs
if opt == 1
for j = 1:10000
    % align time simulations so that CIs post maximum infection start time
    % can be obtained 
    tTest = [Start,t(ind)];
    tt = tTest(1):.1:20;
    tt1 = 0;
    if Starts(j) < Start
    tt1 = linspace(Starts(j),Start-(Start-Starts(j))/2,50);
    tt2 = tTest(1):.1:20;
    tt = [tt1,tt2];
    end
    % get the model solution for a particular parameter set and given host
    p = params(j,:,its);

    
    
        lambda = (1/dec)*min(x)*.8;


    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); % use disperse here ...
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
 
  
   if mod(j,100) == 0
       % print progress of solution generation 
       fprintf('iteration %i of %i\n',j,10000)
   end
end
remove = 250;

% get CIs and median time trajectory from simulations 
% sort the simulations at each time point
sols1 = sort(sols1);
sols2 = sort(sols2);
sols3 = sort(sols3);


% get end points of confidence intervals for plotting 
lb1 = sols1(remove+1,:);
ub1 = sols1(end-remove,:);
lb2 = sols2(remove+1,:);
ub2 = sols2(end-remove,:);
lb3 = sols3(remove+1,:);
ub3 = sols3(end-remove,:);

base(:,1) = median(sols1);
base(:,2) = median(sols2);
base(:,3) = median(sols3);

% strings for generating plot title and save file name respectively 
str = strcat('ID',' ',num2str(IDs(its)),{', '},extractBetween(groups(its),1,4));
str2= strcat('ID',num2str(IDs(its)),groups(its),'CIHapto');


tt = Start:.1:20;
% plot everything 
figure
hold on

% plot model solutions for innate 
x2 = [tt,fliplr(tt)];
inBetween = [lb1,fliplr(ub1)];
h = fill(x2, inBetween,cmap(1,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
scatter(t2,x,120,cmap(1,:),'Marker','+','linewidth',2)

% plot model solutions for virus 
inBetween = [lb2,fliplr(ub2)];
h = fill(x2, inBetween,cmap(3,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(tt,base(:,2),'Color',cmap(3,:),'LineStyle','-','linewidth',2)
scatter(t2,sesC(its,:),120,cmap(3,:),'Marker','*','linewidth',2)

% plot model solutions for adaptive 
inBetween = [lb3,fliplr(ub3)];
h = fill(x2, inBetween,cmap(2,:),'LineStyle','none','marker','none');
set(h,'facealpha',.4)
scatter(t2,pvsC(its,:),120,cmap(2,:),'Marker','o','linewidth',2)
plot(tt,base(:,3),'Color',cmap(2,:),'LineStyle','--','marker','none','linewidth',2)

    fever = readtable('QualCompInd.csv');
    fstart = table2array(fever(its,end-2))-2;
    fend = table2array(fever(its,end-1))-2;
    if fstart > 0
    xf = [fstart,fend];
    xf = [xf,fliplr(xf)];
    inBetween = [[0,0],fliplr([11,11])];
    h = fill(xf, inBetween,'y','LineStyle','none','marker','none');
    set(h,'facealpha',.2)
   
    end
xline(Start,'k--')
% plot simple random sample of solutions prior to max infection start time 
y = randsample(1:10000,20);
for j = 1:length(y)
    tTest = [Start,t(ind)];
    tt = tTest(1):.1:15;
    tt1 = 0;
    if Starts(y(j)) < Start
    tt1 = linspace(Starts(y(j)),Start,50);
       p = params(y(j),:,its);

    
    
        lambda = (1/dec)*min(x).*.8;

 
    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); % use disperse here ...
    x0 = [p(7),p(9),p(8)]; % initial conditions 
    atilde=0;

     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,tt1,x0,options);

    end
    plot(tt1,yy(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',2)
    plot(tt1,yy(:,2),'Color',cmap(3,:),'LineStyle',':','linewidth',2)
    plot(tt1,yy(:,3),'Color',cmap(2,:),'LineStyle',':','linewidth',2)
end
hold off

% plot formatting 
ylabel('Concentration')
ylim([0 11])

 xlim([0 20])

xlabel('Time since contact was initiated (days)')
title(str)
set(gca,'FontSize',16)
% save the plot
baseFileName = sprintf(str2);
% change the string as appropriate for you computer 
fname = '~/Documents/MATLAB/Model_identifiability/Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');


elseif opt == 2 % alternative plotting option, 
    % useful for checking coverage of the empirical data 

% colormaps for differnt serotypes 
cmap1 = [.5 .5 .5 .2;
0 0 0 .2;
0 0 1 .4];
cmap2 = [.5 .5 .5 .2;
0 0 0 .2;
0 1 0 .4];
cmap3 = [.5 .5 .5 .2;
0 0 0 .2;
1 0 0 .4];

if its <= 4
    cmap = cmap1;
end
if its >4 && its <=8
    cmap = cmap2;
end
if its >8
    cmap = cmap3;
end

figure
hold on
% random sample of 500 paramter sets 
y = randsample(1:10000,500);
for j = 1:length(y)

p = params(y(j),:,its);
        lambda = (1/dec)*min(x).*.8;

 % solve the ode   
    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); 
    x0 = [p(7),p(9),p(8)]; % initial conditions 
    atilde=0;

     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
t_plot = [Starts(y(j)),t(t>Starts(y(j)))];
[~,yy_plot] = ode45(f,t_plot,x0,options);
% plot the solution for innate
    plot(t_plot,yy_plot(:,1),'Color',cmap(1,:),'LineStyle',':','linewidth',.5)

end

for j = 1:length(y)
% plot the solution for virus 
p = params(y(j),:,its);
        lambda = (1/dec)*min(x).*.8;


    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); 
    x0 = [p(7),p(9),p(8)]; 
    atilde=0;

     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
t_plot = [Starts(y(j)),t(t>Starts(y(j)))];
[~,yy_plot] = ode45(f,t_plot,x0,options);
    plot(t_plot,yy_plot(:,2),'Color',cmap(3,:),'LineStyle',':','linewidth',.5)

end

for j = 1:length(y)
% plot the solution for adaptive 
p = params(y(j),:,its);
        lambda = (1/dec)*min(x).*.8;


    k=p(1); d=1/dec; r=p(5); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(6); 
    x0 = [p(7),p(9),p(8)]; 
    atilde=0;

     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
t_plot = [Starts(y(j)),t(t>Starts(y(j)))];
[~,yy_plot] = ode45(f,t_plot,x0,options);
    plot(t_plot,yy_plot(:,3),'Color',cmap(2,:),'LineStyle',':','linewidth',.5)

end
% plot the data 
plot(t,x,'Color','black','LineStyle','-','Marker','o','linewidth',2)
plot(t,s,'Color','blue','LineStyle','-','Marker','o','linewidth',2)
plot(t,pv,'Color','black','LineStyle','-','Marker','o','linewidth',2)
hold off
end
end