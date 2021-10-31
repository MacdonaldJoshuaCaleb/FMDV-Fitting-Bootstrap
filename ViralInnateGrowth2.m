function [] = ViralInnateGrowth2(index,opt)
close all
ses = readtable('ViremiaContact.csv');
ses = table2array(ses(index,3:end));

xs1 = readtable('HaptoDataContact.csv');
% IDs = xs(:,2);
% IDs = table2array(IDs)
% groups = xs(:,1);
% groups = table2array(groups);
% groups = string(groups);
xs1 = table2array(xs1(index,3:end));
pvsC = readtable('VNTContactData.csv');
pvsC = table2array(pvsC(index,3:end));

    xs2 = readtable('SAADataContact.csv');
    IDs = xs2(:,2);
IDs = table2array(IDs);
groups = xs2(:,1);
groups = table2array(groups);
groups = string(groups);
    xs2 = table2array(xs2(index,3:end));

if opt == 0
    xs = xs2;
end
if opt == 1
    xs = xs1;
end

tt = [0,2, 4, 6, 9,12,28];


virus = ses;
response = xs;
    for j = 1:length(virus)
        if virus(j) > 0
            I = j;
            break
        end
    end
    virus = virus(1:I);
    response = response(1:I);
    tts = tt(1:I);
    if I > 1
        virus = [virus(1),virus(I)];
        response = [response(1),response(I)];
        tts = [tt(1),tt(I)];
    end
    if opt == 0
        d = 1/11.8;
    end
    if opt == 1
        d = 1/21.23;
    end
    
    function dudt = TheModel(t,u,param)
        k = param(1);
        r = param(2);
        theta = param(3);
        
        dudt = [k*u(1)*u(2) + d*(min(response)-u(1)); r*u(2) - theta*u(1)*u(2)];
    end

wt1 = 1/max(response);
wt2 = 1/max(virus);
  function z = objFxn(param,tt)
        %P0 = p(1);
        I0 = response(1);
        P0 = param(4);
        u0 = [I0,P0];
        [t,Sol] = ode45(@(t,u)TheModel(t,u,param),tt,u0);
            z1 = [Sol(1,1),Sol(end,1)];
            z2 = [Sol(1,2),Sol(end,1)];
            z = [z1*wt1,z2*wt2];
  end
vv = ViralGrowth(ses);
lb = [9e-5, vv(2),9e-9, vv(1)];
ub = [12,12,.01,1];
p0 = [6,1.1*vv(2),1e-5,1.1*vv(1)];

[pfit resnorm] = lsqcurvefit(@objFxn,p0,tts,[response*wt1,virus*wt2],lb,ub);
array2table(pfit)
end
