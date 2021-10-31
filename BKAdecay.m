function [pfit] = BKAdecay(index)
xs = readtable('BKAPlasmaContact.csv');
xs = table2array(xs(:,4:end));
ses = readtable('ViremiaContact.csv');
ses = table2array(ses(:,3:end));
s = ses(index,1:end-1);
%xEcoli = xs(index,1:end-1);
x = xs(index+12,1:end-2);
 if index == 1
     x = x(1:end-1);
 end
for j = 1:length(x)
  %  if xEcoli(j) < 0
   %     xEcoli(j) = 0;
   % end
    if x(j) < 0
        x(j) = 0;
    end
end
tt = [0,2, 4, 6, 9,12,28];


function yy = paramfun1(p,t)
     k = p(1);
     x0 = p(2);
 f = @(t,u) -k*u;
 [~,yy] = ode45(f,t,x0);
 yy = yy';
 if length(tt(ind:stop)) <= 2
     yy = [yy(1),yy(end)];
 end
 end 
response = x;
 [I0,ind] = max(response);
 stop = ind+1
 for j = ind+1:length(response)
     if response(j) < response(j-1)
         stop = j;
     end
 end
 p0 = [1/2.5,.5 ];
lb = [0,0];
ub = [1,1];
%paramfun1(p0,tt(ind:stop))
%response(in
 [pfit resnorm] = lsqcurvefit(@paramfun1,p0,tt(ind:stop),response(ind:stop),lb,ub);

%figure 
%hold on
%plot(tt(ind:stop),paramfun1(pfit,tt(ind:stop)))
%plot(tt(ind:stop),response(ind:stop),'b*')
%hold off
end

