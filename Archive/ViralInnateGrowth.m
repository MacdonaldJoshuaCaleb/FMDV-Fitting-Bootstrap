function [] = ViralInnateGrowth(index,opt)
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
% figure
% hold on
% yyaxis left
% plot(tt(1:I),virus(1:I),'r*')
% yyaxis right
% plot(tt(1:I),response(1:I),'b*')

         if opt == 0 %'SAA'
             d = 1./11.18; % mean decay rate of SAA protein 
         end
        if opt == 1 %'Hapto'
            d = 1./21.23; % mean decay rate of Hapto protein 
        end 
     v = ViralGrowth(ses);
    function dudt = paramfun1(t,u,p)
        %P0 = p(1);
        r = p(2);
        %I0 = p(3);
        k = p(1);
        theta =0;
%         if opt == 1
%         k0 = p(5);
%         end
        %I0 = p(1);
        %k1 = p(6);
        %u0 = [I0,P0];
       
        %d = 0;
        if opt == 1
        dudt = [k*u(1)*u(2)+d*(min(response)-u(1));... % innate response, formm of activation seems better when using indirect markers of innate response 
                        r*u(2)-(theta)*u(1)*u(2)]; % pathogen 
        end
        if opt == 0
              dudt = [k*u(1)*u(2)+d*(min(response)-u(1));... % innate response, formm of activation seems better when using indirect markers of innate response 
                        r*u(2)-(theta)*u(1)*u(2)]; % pathogen 
        end
    end 
%   
%   
% v
%         
% 
%  virus = ses;
%  response = xs;
%  if opt == 1
%  Ip = diff(xs1(1:end-1)./max(xs1(1:end-1)))./2;
%  end
%  if opt == 0
%      Ip = diff(xs2(1:end-1)./max(xs1(1:end-1)))./2;
%  end
% %  figure
% %  plot(tt(2:end-1),Ip)
% %  
% %     for j = 1:length(virus)
% %         if virus(j) > 0 
% %             I = j+1;
% %             break
% %         end
% %     end
% for j = 1:length(Ip)
%     if Ip > .1
%         I = j+1;
%         break
%     end
% end
%     virus = virus(1:I);
%     response = response(1:I);
%     tts = tt(1:I);
%     if opt == 1
     wt1 = 1/max(xs(1:I));
     wt2 = 1/max(ses(1:I));
%     end
%     if opt == 0
%         wt1 = 1/max(xs(1:I-1));
%         wt2 = 1/max(ses(1:I-1));
%     end
%     
response(1)
    function z = objFxn(p,tt)
        %P0 = p(1);
        I0 = response(1);
        P0 = p(end);
        u0 = [I0,P0];
        [t,Sol] = ode45(@(t,u)paramfun1(t,u,p),tt,u0);
        if I <= 2
            if opt == 1
            z1 = [Sol(1,1),Sol(end,1)];
            z2 = [Sol(1,2),Sol(end-1,1)];
            z = [z1*wt1,z2*wt2];
            end
            if opt == 0
                z1 = [Sol(1,1),Sol(end,1)];
                 z2 = [Sol(1,2),Sol(end,2)];
            end
        end
        if I > 2
            if opt == 1
            z = [Sol(1:end,1)'*wt1,Sol(end,2)'*wt2];
            end
            if opt == 0
                z = [Sol(1:end,1)'*wt1,Sol(end,2)'*wt2];
            end
        end
     end  
% %     tp = linspace(tts(1),tts(end));
% %     % [P0, r, I0, k, theta]
% if response(1) > 0
%      M = min(response(1:end-1));
% end
% if response(1) == 0
%     M = min(response(2:end-1));
% end
%  v = ViralGrowth(ses);
% %d
% 
   lb = [ 0,    v(2)+9e-7,  v(1)];
      ub = [4, 10*v(2)  , 1];
 
%     lb = lb(1:end-1);
%     ub = ub(1:end-1);
     lb(1) = 0;
     ub(1) = 5;
%    %  v(2)
% 
  p0 = [d+.05,mean([lb(3),ub(3)]),1.5*v(1)];
% 
%      p0 = p0(1:end-1);
% 
%     % set rate as d/
%   
% 
%  
%     objFxn(p0,tts)
%     if opt == 1
     [pfit resnorm] = lsqcurvefit(@objFxn,p0,tt(1:I),[response*wt1,virus(end)*wt2],lb,ub);
     array2table(pfit)
      I0 = response(1);
        P0 = pfit(end);
        u0 = [I0,P0];
        [t,Sol] = ode45(@(t,u)paramfun1(t,u,pfit),tt(1:I),u0);
        figure
        yyaxis left
        hold on
        plot(t,Sol(:,1))
        plot(tt(1:I),response(1:I),'b*')
        hold off
        yyaxis right
        hold on
        plot(t,Sol(:,2))
        plot(tt(1:I),virus(1:I),'r*')
        hold off
%     end
%     if opt == 0
%         [pfit resnorm] = lsqcurvefit(@objFxn,p0,tts,[response(1:end-1)*wt1,virus(1:end-1)*wt2],lb,ub);
%     end
%     pfit = [pfit,v(1)];
%     array2table(pfit)
%     ttp2 = linspace(tts(1),tts(end));
%     ttp1 = linspace(tts(1),tts(end-1));
%     if opt == 1
%         u0 = [pfit(1),v(1)];
%         [tp1,Solp1] = ode45(@(t,u)paramfun1(t,u,pfit),ttp1,u0);
%         [tp2,Solp2] = ode45(@(t,u)paramfun1(t,u,pfit),ttp2,u0);
%       str = strcat('ID',' ',num2str(IDs(index)),', ',extractBetween(groups(index),1,4));   
%     red = [1 0 0];
%     blue = [0 0 1];
%     left_color = red;
%     right_color = blue;
% 
%         earlyfig = figure;
%         set(earlyfig,'defaultAxesColorOrder',[left_color; right_color]);
%           yyaxis left
%        hold on
%        plot(tp1,Solp1(:,2),'r','linewidth',2)
%        plot(tt(1:end-1),ses(1:end-1),'rd','MarkerFaceColor','r')
%        plot(tt(1:end-1),pvsC(1:end-1),'gd','MarkerFaceColor','g')
%        hold off
%        ylabel('viral load')
%        ylim([0 10])
%         yyaxis right
%         hold on
%         plot(tp2,Solp2(:,1),'b','linewidth',2)
%         plot(tt(1:end-1),xs(1:end-1),'bd','MarkerFaceColor','b')
%         if opt == 1
%         ylabel('\mug/mL Haptoglobin')
%         end
%         if opt == 0
%             ylabel('\mug/mL SAA')
%         end
%         hold off
%         title(str)
%         xlim([tt(1) tt(end-1)])
%         xlabel('infection age')
%          set(gca,'FontSize',16)
%     end
%       if opt == 0
%         u0 = [pfit(1),v(1)];
%         [tp1,Solp1] = ode45(@(t,u)paramfun1(t,u,pfit),ttp1,u0);
%       %  [tp2,Solp2] = ode45(@(t,u)paramfun1(t,u,pfit),ttp2,u0);
%       str = strcat('ID',' ',num2str(IDs(index)),', ',extractBetween(groups(index),1,4));   
%     red = [1 0 0];
%     blue = [0 0 1];
%     left_color = red;
%     right_color = blue;
% 
%         earlyfig = figure;
%         set(earlyfig,'defaultAxesColorOrder',[left_color; right_color]);
%           yyaxis left
%        hold on
%        plot(tp1,Solp1(:,2),'r','linewidth',2)
%        plot(tt(1:end-1),ses(1:end-1),'rd','MarkerFaceColor','r')
%        plot(tt(1:end-1),pvsC(1:end-1),'gd','MarkerFaceColor','g')
%        hold off
%        ylabel('viral load')
%          ylim([0 10])
%         yyaxis right
%         hold on
%         plot(tp1,Solp1(:,1),'b','linewidth',2)
%         plot(tt(1:end-1),xs(1:end-1),'bd','MarkerFaceColor','b')
%         if opt == 1
%         ylabel('\mug/mL Haptoglobin')
%         end
%         if opt == 0
%             ylabel('\mug/mL SAA')
%         end
%         hold off
%         title(str)
%         xlim([tt(1) tt(end-1)])
%         xlabel('infection age')
%          set(gca,'FontSize',16)
%      end
% %     %[vpfit resnorm] = lsqcurvefit(@paramfun1,vpfit,tts,virus,lb,ub);
%      fit = paramfun1(pfit,tp);
%      Ifit = fit(1:I);
%      Pfit = fit(I+1:end);
%      figure
%      hold on
%      plot(tp,Pfit,'r',tts,virus,'d','color','red')
%      plot(tp,Ifit,'b',tts,response,'d','color','blue')
%      hold off
% %     Virusfits(j,:) = vpfit; 
% %     save('Virusfits.mat','Virusfits')
% % xx = (xs1+xs2)./2;
% % figure
% % plot(tt,xx)
% % title(str)
% 
% % Data = readtable('SortedCleanedBKAContact.csv');
% % EcoliData = table2array(Data(1:12,3:end));
% % StaphData = table2array(Data(13:24,3:end));
% % BKA1 = EcoliData(index,:);
% % BKA2 = StaphData(index,:);
% % for w = 1:length(BKA1)
% %     if BKA1(w) < 0
% %         BKA1(w) = 0;
% %     end
% %     if BKA2(w) < 0
% %         BKA2(w) = 0;
% %     end
% % end
% % figure
% % 
% % hold on
% % plot(tt(1:end-1),xs1(1:end-1)./max(xs1(1:end-1)),'b--d','MarkerFaceColor','b')
% % plot(tt(1:end-1),xs2(1:end-1)./max(xs2(1:end-1)),'b-.^', 'MarkerFaceColor','b')
% % plot(tt(1:end-2),BKA1(1:end-2)./max(BKA1(1:end-2)),'k--d','MarkerFaceColor','k')
% % plot(tt(1:end-2),BKA2(1:end-2)./max(BKA2(1:end-2)),'k:^','MarkerFaceColor','k')
% % xline(tt(I-1),'r');
% % hold off
% % title(str)
% % ylabel('prop. maximum value')
% % legend({'Hapto','SAA','BKA e. coli','BKA staph', 'P(\tau) > 0'},'Location','NorthEastOutside')
% % ylim([0 inf])
% % xlim([0, tt(I+1)])
% % figure
% % plot(tt(2:end-1),diff(xs1(1:end-1))./2)
end