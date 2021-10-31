close all
Data = readtable('SortedCleanedBKAContact.csv');
EcoliData = table2array(Data(1:12,3:end));
StaphData = table2array(Data(13:24,3:end));
IDs = table2array(Data(1:12,2));
ts = [0, 2, 4, 6, 9, 12, 28]; 
colors = jet(12);
ses = readtable('ViremiaContact.csv');
ses = table2array(ses(:,3:end));






for j = 1:12
    figure
   index = find(EcoliData(j,1:end-1) <= 0);
   EcoliData(j,index) = 0;
   str(j) = strcat("ID ",num2str(IDs(j)));
plot(ts(1:end-2),EcoliData(j,1:end-2),'-.*','Color',colors(j,:),'linewidth',2)
end
%legend(str,'NumColumns',3,'location','NorthWest','fontsize',12)
xlim([ts(1), ts(end-2)])
ylim([0,.45])
ylabel('P BKA Ecoli')
hold off
title('Accute Phase')
set(gca,'FontSize',16)
% 
% figure
% hold on
% for j = 1:12
%    index = find(StaphData(j,:) <= 0);
%   StaphData(j,index) = 0;
%    str(j) = strcat("ID ",num2str(IDs(j)));
% plot(ts,StaphData(j,:),'-.*','Color',colors(j,:),'linewidth',2)
% end
% legend(str,'NumColumns',3,'location','NorthWest','fontsize',12)
% xlim([ts(1), ts(end)])
% ylim([0, 1.5])
% ylabel('P BKA Staph')
% hold off
% title('Full Timeframe')
% set(gca,'FontSize',16)



% for j = 1:12
%     figure
%    index = find(StaphData(j,1:end-1) <= 0);
%    StaphData(j,index) = 0;
%    str(j) = strcat("ID ",num2str(IDs(j)));
% plot(ts(1:end-1),StaphData(j,1:end-1),'-.*','Color',colors(j,:),'linewidth',2)
% end
% legend(str,'location','NorthWest','fontsize',12)
% xlim([ts(1), ts(end-1)])
% %ylim([0,1.5])
% ylabel('P BKA Staph')
% hold off
% title('Accute Phase')
% set(gca,'FontSize',16)