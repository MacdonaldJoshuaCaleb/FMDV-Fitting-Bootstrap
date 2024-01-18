function pfit = IFNaDecay(xs)
%ses = readtable('ViremiaContact.csv');
%ses = table2array(ses(index,3:end))
tt = [0,2, 4, 6, 9,12,28];




 function yy = paramfun1(p,t)
     k = p(1);
     x0 = p(2);
 f = @(t,u) -k*u;
 [~,yy] = ode45(f,t,x0);
 yy = yy';
 if length(tt(ind:end)) <= 2
     yy = [yy(1),yy(end)];
 end
 end 
response = xs;
 [I0,ind] = max(response);
 p0 = [1/2.5, I0];
lb = [0,.8*I0];
ub = [1,1.2*I0];
paramfun1(p0,tt(ind:end))
response(ind:end)
 [pfit resnorm] = lsqcurvefit(@paramfun1,p0,tt(ind:end),response(ind:end),lb,ub);

%figure 
%plot(tt(ind,end))
%hold on
%plot(tt(ind:end),paramfun1(pfit,tt(ind:end)))
%plot(tt(ind:end),response(ind:end),'b*')
%hold off
end