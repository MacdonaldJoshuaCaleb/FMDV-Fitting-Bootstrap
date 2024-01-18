function [] = GenerateCI
close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsC = readtable('pvfitsC.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsC.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsC.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('SAADataContact.csv');
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
pvsN = readtable('VNTNeedleData.csv');function [] = GenerateCI
close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsC = readtable('pvfitsC.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsC.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsC.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('SAADataContact.csv');
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

Fits1 = zeros(its,3*length(tt));
Fits2 = zeros(its,3*length(tt));
Fits3 = zeros(its,3*length(tt));
Fits4 = zeros(its,3*length(tt));

params1 = zeros(its,7);
params2 = zeros(its,7);
params3 = zeros(its,7);
params4 = zeros(its,7);

cumviralContact1 = zeros(its,1);
cumviralContact2 = zeros(its,1);
cumviralContact3 = zeros(its,1);
cumviralContact4 = zeros(its,1);

maxViralContact1 = zeros(its,2);
maxViralContact2 = zeros(its,2);
maxViralContact3 = zeros(its,2);
maxViralContact4 = zeros(its,2);

SData1 = zeros(its,21);
SData2 = zeros(its,21);
SData3 = zeros(its,21);
SData4 = zeros(its,21);

function yy = paramfun1(p,t)
         lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(:,1)',yy(:,2)',yy(:,3)'];
end 

for count = 1:its
for index = 1:4
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 1
    params1(count,:) = FittingforCI(dd,index);
    Fits1(count,:) = paramfun1(params1(count,:),tt);
    [mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    maxViralContact1(count,:) = [mVi tt(I)];
    SData1(count,:) = dd;
end
if index == 2
    params2(count,:) = FittingforCI(dd,index);
    Fits2(count,:) = paramfun1(params2(count,:),tt);function [] = GenerateCI
close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
pvfitsC = readtable('pvfitsC.csv');
pvfitsC = table2array(pvfitsC);
%pvfitsN = readtable('pvfitsN.csv');
%pvfitsN = table2array(pvfitsN);
sfitsC = readtable('sfitsC.csv');
sfitsC = table2array(sfitsC);
%sfitsN = readtable('sfitsN.csv');
%sfitsN = table2array(sfitsN);
xfitsC = readtable('xfitsC.csv');
xfitsC = table2array(xfitsC);
%xfitsN = readtable('xfitsN.csv');
%xfitsN = table2array(xfitsN);
xsC = readtable('SAADataContact.csv');
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

Fits1 = zeros(its,3*length(tt));
Fits2 = zeros(its,3*length(tt));
Fits3 = zeros(its,3*length(tt));
Fits4 = zeros(its,3*length(tt));

params1 = zeros(its,7);
params2 = zeros(its,7);
params3 = zeros(its,7);
params4 = zeros(its,7);

cumviralContact1 = zeros(its,1);
cumviralContact2 = zeros(its,1);
cumviralContact3 = zeros(its,1);
cumviralContact4 = zeros(its,1);

maxViralContact1 = zeros(its,2);
maxViralContact2 = zeros(its,2);
maxViralContact3 = zeros(its,2);
maxViralContact4 = zeros(its,2);

SData1 = zeros(its,21);
SData2 = zeros(its,21);
SData3 = zeros(its,21);
SData4 = zeros(its,21);

function yy = paramfun1(p,t)
         lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(:,1)',yy(:,2)',yy(:,3)'];
end 

for count = 1:its
for index = 1:4
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 1
    params1(count,:) = FittingforCI(dd,index);
    Fits1(count,:) = paramfun1(params1(count,:),tt);
    [mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    maxViralContact1(count,:) = [mVi tt(I)];
    SData1(count,:) = dd;
end
if index == 2
    params2(count,:) = FittingforCI(dd,index);
    Fits2(count,:) = paramfun1(params2(count,:),tt);
    [mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    maxViralContact2(count,:) = [mVi tt(I)];
    SData2(count,:) = dd;
end
if index == 3
    params3(count,:) = FittingforCI(dd,index);
    Fits3(count,:) = paramfun1(params3(count,:),tt);
    [mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    maxViralContact3(count,:) = [mVi tt(I)];
    SData3(count,:) = dd;
end
if index == 4
    params4(count,:) = FittingforCI(dd,index);
    Fits4(count,:) = paramfun1(params4(count,:),tt);
    [mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    maxViralContact4(count,:) = [mVi tt(I)];
    SData4(count,:) = dd;
end

end
fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end



save('Fits1.mat','Fits1')
save('Fits2.mat','Fits2')
save('Fits3.mat','Fits3')
save('Fits4.mat','Fits4')

save('params1.mat','params1')
save('params2.mat','params2')
save('params3.mat','params3')
save('params4.mat','params4')

save('cumviralContact1.mat','cumviralContact1')
save('cumviralContact2.mat','cumviralContact2')
save('cumviralContact3.mat','cumviralContact3')
save('cumviralContact4.mat','cumviralContact4')

save('maxViralContact1.mat','maxViralContact1')
save('maxViralContact2.mat','maxViralContact2')
save('maxViralContact3.mat','maxViralContact3')
save('maxViralContact4.mat','maxViralContact4')

save('SData1.mat','SData1')
save('SData2.mat','SData2')
save('SData3.mat','SData3')
save('SData4.mat','SData4')
end
    [mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    maxViralContact2(count,:) = [mVi tt(I)];
    SData2(count,:) = dd;
end
if index == 3
    params3(count,:) = FittingforCI(dd,index);
    Fits3(count,:) = paramfun1(params3(count,:),tt);
    [mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    maxViralContact3(count,:) = [mVi tt(I)];
    SData3(count,:) = dd;
end
if index == 4
    params4(count,:) = FittingforCI(dd,index);
    Fits4(count,:) = paramfun1(params4(count,:),tt);
    [mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    maxViralContact4(count,:) = [mVi tt(I)];
    SData4(count,:) = dd;
end

end
fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end



save('Fits1.mat','Fits1')
save('Fits2.mat','Fits2')
save('Fits3.mat','Fits3')
save('Fits4.mat','Fits4')

save('params1.mat','params1')
save('params2.mat','params2')
save('params3.mat','params3')
save('params4.mat','params4')

save('cumviralContact1.mat','cumviralContact1')
save('cumviralContact2.mat','cumviralContact2')
save('cumviralContact3.mat','cumviralContact3')
save('cumviralContact4.mat','cumviralContact4')

save('maxViralContact1.mat','maxViralContact1')
save('maxViralContact2.mat','maxViralContact2')
save('maxViralContact3.mat','maxViralContact3')
save('maxViralContact4.mat','maxViralContact4')

save('SData1.mat','SData1')
save('SData2.mat','SData2')
save('SData3.mat','SData3')
save('SData4.mat','SData4')
end
pvsN = table2array(pvsN(:,3:end));

tt = 0:.1:30;
its = 10000;

Fits1 = zeros(its,3*length(tt));
Fits2 = zeros(its,3*length(tt));
Fits3 = zeros(its,3*length(tt));
Fits4 = zeros(its,3*length(tt));

params1 = zeros(its,7);
params2 = zeros(its,7);
params3 = zeros(its,7);
params4 = zeros(its,7);

cumviralContact1 = zeros(its,1);
cumviralContact2 = zeros(its,1);
cumviralContact3 = zeros(its,1);
cumviralContact4 = zeros(its,1);

maxViralContact1 = zeros(its,2);
maxViralContact2 = zeros(its,2);
maxViralContact3 = zeros(its,2);
maxViralContact4 = zeros(its,2);

SData1 = zeros(its,21);
SData2 = zeros(its,21);
SData3 = zeros(its,21);
SData4 = zeros(its,21);

function yy = paramfun1(p,t)
         lambda=0; k=p(3); d=1/11.18; r=p(1); K=p(4); theta=0; delta=p(5); b=p(6); % use disperse here ...
    x0 = [p(7),p(2),.0001*p(7)]; % initial conditions 
    atilde=0;
    f = @(t,a) [lambda+k*a(2)*a(1)- d*a(1); ... 
                (r*(1-(a(2)/K))-theta*a(1)-delta*a(3))*a(2); ...
                (atilde*(a(1)/(1+a(1)))+b*a(3))*a(2)];
   
   [~,yy] = ode45(f,t,x0,options);
   yy = [yy(:,1)',yy(:,2)',yy(:,3)'];
end 

for count = 1:its
for index = 1:4
nLevel = [.4,.4,.4];
dd = GetData(xfitsC(index,:),sfitsC(index,:),pvfitsC(index,:),nLevel);

if index <= 12    % contact
    t = [0,2, 4, 6, 9, 12, 28];
end
if index > 12    % needle
    t = [0,2, 4, 6, 8, 11, 14, 30];
end
if index == 1
    params1(count,:) = FittingforCI(dd,index);
    Fits1(count,:) = paramfun1(params1(count,:),tt);
    [mVi I] = max(Fits1(length(tt)+1:2*length(tt)));
    cumviralContact1(count) = sum(Fits1(count,length(tt)+1:2*length(tt)));
    maxViralContact1(count,:) = [mVi tt(I)];
    SData1(count,:) = dd;
end
if index == 2
    params2(count,:) = FittingforCI(dd,index);
    Fits2(count,:) = paramfun1(params2(count,:),tt);
    [mVi I] = max(Fits2(length(tt)+1:2*length(tt)));
    cumviralContact2(count) = sum(Fits2(count,length(tt)+1:2*length(tt)));
    maxViralContact2(count,:) = [mVi tt(I)];
    SData2(count,:) = dd;
end
if index == 3
    params3(count,:) = FittingforCI(dd,index);
    Fits3(count,:) = paramfun1(params3(count,:),tt);
    [mVi I] = max(Fits3(length(tt)+1:2*length(tt)));
    cumviralContact3(count) = sum(Fits3(count,length(tt)+1:2*length(tt)));
    maxViralContact3(count,:) = [mVi tt(I)];
    SData3(count,:) = dd;
end
if index == 4
    params4(count,:) = FittingforCI(dd,index);
    Fits4(count,:) = paramfun1(params4(count,:),tt);
    [mVi I] = max(Fits4(length(tt)+1:2*length(tt)));
    cumviralContact4(count) = sum(Fits4(count,length(tt)+1:2*length(tt)));
    maxViralContact4(count,:) = [mVi tt(I)];
    SData4(count,:) = dd;
end

end
fprintf('-----------------------\n')
fprintf('Iteration %i of %i\n', 4*count,4*its)
fprintf('-----------------------\n')
end



save('Fits1.mat','Fits1')
save('Fits2.mat','Fits2')
save('Fits3.mat','Fits3')
save('Fits4.mat','Fits4')

save('params1.mat','params1')
save('params2.mat','params2')
save('params3.mat','params3')
save('params4.mat','params4')

save('cumviralContact1.mat','cumviralContact1')
save('cumviralContact2.mat','cumviralContact2')
save('cumviralContact3.mat','cumviralContact3')
save('cumviralContact4.mat','cumviralContact4')

save('maxViralContact1.mat','maxViralContact1')
save('maxViralContact2.mat','maxViralContact2')
save('maxViralContact3.mat','maxViralContact3')
save('maxViralContact4.mat','maxViralContact4')

save('SData1.mat','SData1')
save('SData2.mat','SData2')
save('SData3.mat','SData3')
save('SData4.mat','SData4')
end

% need to bring in time to max viremia, max viremia, and viremia duration 
% can get cum