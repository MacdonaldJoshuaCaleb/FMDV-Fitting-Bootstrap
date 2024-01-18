close all
%clear
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
t1 = [0,2, 4, 6, 9,12,28];

set(0,'DefaultFigureVisible','on')
% import and organize data
xs = readtable('HaptoDataContact.csv');

Starts = readtable('InfectionStartTimes.csv');
Starts = table2array(Starts);
%Starts = mean(Starts);

groups = xs(:,1);
groups = table2array(groups);
groups = string(groups);

IDs = xs(:,2);
IDs = table2array(IDs);
%IDs = string(groups);

xs = table2array(xs(:,3:end));
xs(xs>0) = log10(xs(xs>0));

ses = readtable('ViremiaContact.csv');
ses = table2array(ses(:,3:end));

pvs = readtable('VNTContactData.csv');
pvs = table2array(pvs(:,3:end));


    
cmap1 = [.5 .5 .5;
0 0 0;
0 0 1];
cmap2 = [.5 .5 .5;
0 0 0;
0 1 0];
cmap3 = [.5 .5 .5;
0 0 0;
1 0 0];
   %yy = [yy(1:end-1,1)'.*sqwt3,yy(1:end-1,2)'.*sqwt,yy(:,3)'.*sqwt2];

load('paramsNoiseAdjust.mat')
params = paramsNoise;
its = [1:9,11:12];
%its = 5;
count = 10000;
for kk = 1:length(its)
    fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
    fprintf('WORKING ON HOSt %i of %i\n',kk,length(its))
    fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
if its(kk) <= 4
    cmap = cmap1;
end
if its(kk) >4 && its(kk) <=8
    cmap = cmap2;
end
if its(kk) >8
    cmap = cmap3;
end
for j = 1:count
%temp = RevisedFittingContactDataHaptoStimes(its(kk),1);
%RE(j,:,kk) = 100*temp(3,:);
%params(j,:,kk) = temp(1,:);
%paramsNoise(j,:,kk) = temp(2,:);
if mod(j,count/100) == 0
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('PERCENT COMPLETE: %i\n',round((j/count)*100))
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
% fprintf('ooooooooooooooooooooooooooooooooooooooo\n')
end
end
%ARE(kk,:) = mean(RE(:,:,kk));
x = xs(its(kk),:);
s = ses(its(kk),:);
pv = pvs(its(kk),:);
    str = strcat('ID',' ',num2str(IDs(its(kk))),', ',extractBetween(groups(its(kk)),1,4));
    str2= strcat('ID',num2str(IDs(its(kk))),groups(its(kk)),'ST','Hapto');
% close all
% figure
% hold on
MST = max(params(:,12,kk));


for jj = 1:10000
    if mod(jj,100) == 0
       fprintf('iteration %i of %i\n',jj,10000)
   end
    
    t = linspace(MST,15,100);
    p = params(jj,1:8,kk);

    if params(jj,12,kk) < MST
        tt = linspace(MST,15,100);

       
        t = [linspace(params(jj,12,kk),MST-(MST-params(jj,12,kk))/20,20),tt];
    end
    temp = [linspace(params(jj,12,kk),MST-(MST-params(jj,12,kk))/20,20),MST];
    tEarly(jj,:) = temp;
    
    
    
        lambda1=(1/21.23)*p(end-2);
    
        lambda2 = (1/21.23)*x(end);
        lambda = min(lambda1,lambda2);
    %lambda = lambda2;
    k=p(1); d=1/21.23; r=p(6); K=p(2); theta=0; delta=p(3); b=p(4); v0 = p(7); % use disperse here ...
    x0 = [p(5),p(8),1e-4*p(5)]; % initial conditions 
    atilde=0;
 %   f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
  %              (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
   %             (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
     f = @(t,a) [lambda+k*(a(2)/(v0+a(2)))*a(1)- d*a(1); ... 
               (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
%     plot(t,yy(:,1),'b')
%     plot(t,yy(:,2),'r')
%     plot(t,yy(:,3),'g')
    if params(jj,12,kk) < MST
    Innate(jj,:) = yy(21:end,1)';
    Virus(jj,:) = yy(21:end,2)';
    Adapt(jj,:) = yy(21:end,3)';
    Iearly(jj,:) = yy(1:21,1)';
    Vearly(jj,:) = yy(1:21,2)';
    Aearly(jj,:) = yy(1:21,3)';
    end
    if params(jj,12,kk) == MST
    Innate(jj,:) = yy(1:end,1)';
    Virus(jj,:) = yy(1:end,2)';
    Adapt(jj,:) = yy(1:end,3)';
    end
end
for i = 1:10000
    for m = 1:21
        if Iearly(i,m) == 0
            Iearly(i,m) = NaN;
            Vearly(i,m) = NaN;
            Aearly(i,m) = NaN;
        end
    end
end
Innate = sort(Innate);
Virus = sort(Virus);
Adapt = sort(Adapt);
ILB = Innate(251,:);
IUB = Innate(end-250,:);
VLB = Virus(251,:);
VUB = Virus(end-250,:);
ALB = Adapt(251,:);
AUB = Adapt(end-250,:);

t = linspace(MST,15,100);

x2 = [t,fliplr(t)];
close all
figure
hold on
inBetween = [ILB,fliplr(IUB)];
h = fill(x2, inBetween,cmap(1,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(t,median(Innate),'Color',cmap(1,:),'linewidth',2)
inBetween = [VLB,fliplr(VUB)];
h = fill(x2, inBetween,cmap(3,:),'LineStyle','none');
set(h,'facealpha',.4)
plot(t,median(Virus),'Color',cmap(3,:),'linewidth',2)
inBetween = [ALB,fliplr(AUB)];
h = fill(x2, inBetween,cmap(2,:),'LineStyle','none');
set(h,'facealpha',.4)
rr = randi([1,10000],20,1);
for zz = 1:length(rr)
    plot(tEarly(rr(zz),:),Iearly(rr(zz),:),'LineStyle',':','Color',cmap(1,:),'linewidth',2)
    plot(tEarly(rr(zz),:),Vearly(rr(zz),:),'LineStyle',':','Color',cmap(3,:),'linewidth',2)
    
 
end
for zz = 1:length(rr)
    plot(tEarly(rr(zz),:),Aearly(rr(zz),:),'Color',cmap(2,:),'linewidth',2)
end
plot(t,median(Adapt),'Color',cmap(2,:),'linewidth',2)
scatter(t1,x,240,'k','Marker','s','linewidth',2)
scatter(t1,s,240,cmap(3,:),'Marker','^','linewidth',2)
scatter(t1,pv,240,cmap(2,:),'Marker','o','linewidth',2)
xline(MST,'k--','linewidth',2)
 fever = readtable('QualCompInd.csv');
    fstart = table2array(fever(its(kk),15))-2;
    fend = table2array(fever(its(kk),16))-2;
    if fstart > 0
    xf = [fstart,fend];
    xf = [xf,fliplr(xf)];
    inBetween = [[0,0],fliplr([11,11])];
    h = fill(xf, inBetween,[0.9290 0.6940 0.1250],'LineStyle','none','marker','none');
    set(h,'facealpha',.2)
    %text(fstart,9,"fever detected")
    end
% plot(t1,s,'+','MarkerSize',20, 'linewidth',2)
% plot(t1,pv,'o','MarkerSize',20,'linewidth',2)
ylabel('Concentration')
xlabel('Time since contact was initiated (days)')
xlim([0 15])
ylim([0 11])
 hold off
title(str)
set(gca,'FontSize',16)
baseFileName = sprintf(str2);
        fname = '~/Documents/MATLAB/Model_identifiability/Plots';
    saveas(gca, fullfile(fname, baseFileName), 'png');
fprintf('DONE\n')
%array2table(ARE)
%  save('AREAdjust.mat','ARE')
%  save('paramsFinalNewAdjust.mat','params')
%  save('REAdjust.mat','RE')
%  save('paramsNoiseAdjust.mat','paramsNoise')

end


