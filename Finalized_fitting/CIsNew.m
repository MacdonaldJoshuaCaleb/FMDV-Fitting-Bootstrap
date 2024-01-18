% generate 95% CIs for model parameters - contact infected hosts 
load('params_noise.mat')
digits = 32;
params = params_noise;
[n k] = size(params(:,:,1));
SAT1Means = zeros(n,k);
SAT2Meeans = zeros(n,k);
SAT3Means = zeros(n,k);
SAT1Meds = zeros(n,k);
SAT2Meds = zeros(n,k);
SAT3Meds = zeros(n,k);
DM13 = zeros(n,k);
DM12 = zeros(n,k);
DM23 = zeros(n,k);
DMed13 = zeros(n,k);
DMed12 = zeros(n,k);
DMed23 = zeros(n,k);

for j = 1:length(params(:,1,1))
SAT1Means(j,:) = (params(j,:,1)+params(j,:,2)+params(j,:,3)+params(j,:,4))./4;
SAT2Means(j,:) = (params(j,:,5)+params(j,:,6)+params(j,:,7)+params(j,:,8))./4;
SAT3Means(j,:) = (params(j,:,9)+params(j,:,10)+params(j,:,11)+params(j,:,12))./4;
SAT1Meds(j,:) = median([params(j,:,1);params(j,:,2);params(j,:,3);params(j,:,4)]);
SAT2Meds(j,:) = median([params(j,:,5);params(j,:,6);params(j,:,7);params(j,:,8)]);
SAT3Meds(j,:) = median([params(j,:,9);params(j,:,10);params(j,:,11);params(j,:,12)]);
DMed13(j,:) = SAT1Meds(j,:)-SAT3Meds(j,:);
DMed12(j,:) = SAT1Meds(j,:)-SAT2Meds(j,:);
DMed23(j,:) = SAT2Meds(j,:)-SAT3Meds(j,:);
DM13(j,:) = (SAT1Means(j,:) - SAT3Means(j,:));
DM12(j,:) = (SAT1Means(j,:) - SAT2Means(j,:));
DM23(j,:) = (SAT2Means(j,:) - SAT3Means(j,:));
end
remove = 250;
DM13 = sort(DM13);
DM12 = sort(DM12);
DM23 = sort(DM23);
DMed13 = sort(DMed13);
DMed12 = sort(DMed12);
DMed23 = sort(DMed23);
T = [SAT1Means;SAT2Means;SAT3Means];
SAT1Means = sort(SAT1Means);
SAT2Means = sort(SAT2Means);
SAT3Means = sort(SAT3Means);
fprintf('Diff Means\n')
AA=[SAT1Means(remove+1,:);SAT1Means(end-remove,:)]
BB=[SAT2Means(remove+1,:);SAT2Means(end-remove,:)]
CC=[SAT3Means(remove+1,:);SAT3Means(end-remove,:)]
fprintf('--------------------------------\n')
% fprintf('Diff Medians')
A=[DM13(remove+1,:);DM13(end-remove,:)]
B=[DM12(remove+1,:);DM12(end-remove,:)]
C=[DM23(remove+1,:);DM23(end-remove,:)]

center_13 = mean(DMed13);
center_12 = mean(DMed12);
center_23 = mean(DMed23);

for i = 1:13
    count = 0;
for j = 1:10000
    if center_13(i) > 0
        if DMed13(j,i) > 0
            count = count+1;
        end
    elseif center_13(i) < 0
        if DMed13(j,i) < 0
            count=count+1;
        end
    end
end
p_13(i) = 1-count/10000;
end



for i = 1:13
    count = 0;
for j = 1:10000
    if center_12(i) > 0
        if DMed12(j,i) > 0
            count = count+1;
        end
    elseif center_12(i) < 0
        if DMed12(j,i) < 0
            count=count+1;
        end
    end
end
p_12(i) = 1-count/10000;
end

for i = 1:13
    count = 0;
for j = 1:10000
    if center_23(i) > 0
        if DMed23(j,i) > 0
            count = count+1;
        end
    elseif center_23(i) < 0
        if DMed23(j,i) < 0
            count=count+1;
        end
    end
end
p_23(i) = 1-count/10000;
end
    
        

digits = 32;
mns(1,:) = mean(DM12(:,1:13));
mns(2,:) = mean(DM13(:,1:13));
mns(3,:) = mean(DM23(:,1:13));
lbs = [B(1,:)',A(1,:)',C(1,:)'];
ubs = [B(2,:)',A(2,:)',C(2,:)'];
s_table_diff= strings([13,3]);
for i = 1:size(s_table_diff,1)
    for j = 1:size(s_table_diff,2)
        s_table_diff(i,j) = strcat(num2str(mns(j,i),'%.3g'),{', ('},num2str(lbs(i,j),'%.3g'),{', '},num2str(ubs(i,j),'%.3g'),')');
    end
end

s_table_diff

mns(1,:) = mean(SAT1Means(:,1:13));
mns(2,:) = mean(SAT2Means(:,1:13));
mns(3,:) = mean(SAT3Means(:,1:13));
lbs = [AA(1,:)',BB(1,:)',CC(1,:)'];
ubs = [AA(2,:)',BB(2,:)',CC(2,:)'];
s_table= strings([13,3]);
for i = 1:size(s_table_diff,1)
    for j = 1:size(s_table_diff,2)
        s_table(i,j) = strcat(num2str(mns(j,i),'%.3g'),{', ('},num2str(lbs(i,j),'%.3g'),{', '},num2str(ubs(i,j),'%.3g'),')');
    end
end


s_table
