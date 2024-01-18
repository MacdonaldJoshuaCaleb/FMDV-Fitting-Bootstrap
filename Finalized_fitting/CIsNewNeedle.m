load('params_noise_needle.mat')
digits = 32;
params = params_noise_needle;
[n k] = size(params(:,:,1));
SAT1Means_needle = zeros(n,k);
SAT2Means_needle = zeros(n,k);
SAT3Means_needle = zeros(n,k);
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
SAT1Means_needle(j,:) = (params(j,:,1)+params(j,:,2)+params(j,:,3)+params(j,:,4))./4;
SAT2Means_needle(j,:) = (params(j,:,5)+params(j,:,6)+params(j,:,7)+params(j,:,8))./4;
SAT3Means_needle(j,:) = (params(j,:,9)+params(j,:,10)+params(j,:,11)+params(j,:,12))./4;
SAT1Meds(j,:) = median([params(j,:,1);params(j,:,2);params(j,:,3);params(j,:,4)]);
SAT2Meds(j,:) = median([params(j,:,5);params(j,:,6);params(j,:,7);params(j,:,8)]);
SAT3Meds(j,:) = median([params(j,:,9);params(j,:,10);params(j,:,11);params(j,:,12)]);
DMed13(j,:) = SAT1Meds(j,:)-SAT3Meds(j,:);
DMed12(j,:) = SAT1Meds(j,:)-SAT2Meds(j,:);
DMed23(j,:) = SAT2Meds(j,:)-SAT3Meds(j,:);
DM13(j,:) = (SAT1Means_needle(j,:) - SAT3Means_needle(j,:));
DM12(j,:) = (SAT1Means_needle(j,:) - SAT2Means_needle(j,:));
DM23(j,:) = (SAT2Means_needle(j,:) - SAT3Means_needle(j,:));
end
remove = 250;
DM13 = sort(DM13);
DM12 = sort(DM12);
DM23 = sort(DM23);
DMed13 = sort(DMed13);
DMed12 = sort(DMed12);
DMed23 = sort(DMed23);
T = [SAT1Means_needle;SAT2Means_needle;SAT3Means_needle];
SAT1Means_needle = sort(SAT1Means_needle);
SAT2Means_needle = sort(SAT2Means_needle);
SAT3Means_needle = sort(SAT3Means_needle);
fprintf('Diff Means\n')
AA=[SAT1Means_needle(remove+1,:);SAT1Means_needle(end-remove,:)]
BB=[SAT2Means_needle(remove+1,:);SAT2Means_needle(end-remove,:)]
CC=[SAT3Means_needle(remove+1,:);SAT3Means_needle(end-remove,:)]
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

load('params_noise.mat')
digits = 32;
params = params_noise;
[n k] = size(params(:,:,1));
SAT1Means = zeros(n,k);
SAT2Means = zeros(n,k);
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
        
% 
% for j = 1:11
%     temp = sort(params(:,:,j));
%     fprintf('CIs\n')
%     lbs(j) = temp(remove+1,:);
%     ubs(j) = temp(end-remove,:);
%     fprintf('Mean\n')
%     mean(temp)
%     fprintf('Median\n')
%     median(temp)
% end

% for j = 1:11
%     temp2(j,:) = median(params(:,:,j));
%     fprintf('CIs\n')
%     array2table([temp(remove+1,:);temp(end-remove,:)])
%     fprintf('Mean\n')
%     mean(temp)
%     fprintf('Median\n')
%     median(temp)
% end
digits = 32;
mns(1,:) = mean(DM12(:,1:13));
mns(2,:) = mean(DM13(:,1:13));
mns(3,:) = mean(DM23(:,1:13));
lbs = [B(1,:)',A(1,:)',C(1,:)'];
ubs = [B(2,:)',A(2,:)',C(2,:)'];
s_table_diff= strings([12,3]);
for i = 1:size(s_table_diff,1)
    for j = 1:size(s_table_diff,2)
        s_table_diff(i,j) = strcat(num2str(mns(j,i),'%.3g'),{', ('},num2str(lbs(i,j),'%.3g'),{', '},num2str(ubs(i,j),'%.3g'),')');
    end
end

s_table_diff

mns(1,:) = mean(SAT1Means_needle(:,1:13));
mns(2,:) = mean(SAT2Means_needle(:,1:13));
mns(3,:) = mean(SAT3Means_needle(:,1:13));
lbs = [AA(1,:)',BB(1,:)',CC(1,:)'];
ubs = [AA(2,:)',BB(2,:)',CC(2,:)'];
s_table= strings([12,3]);
for i = 1:size(s_table_diff,1)
    for j = 1:size(s_table_diff,2)
        s_table(i,j) = strcat(num2str(mns(j,i),'%.3g'),{', ('},num2str(lbs(i,j),'%.3g'),{', '},num2str(ubs(i,j),'%.3g'),')');
    end
end


s_table

DCN1 = SAT1Means-SAT1Means_needle;
DCN2 = SAT2Means-SAT2Means_needle;
DCN3 = SAT3Means-SAT3Means_needle;

DCN1 = sort(DCN1);
DCN2 = sort(DCN2);
DCN3 = sort(DCN3);

A=[DCN1(remove+1,:);DCN1(end-remove,:)]
B=[DCN2(remove+1,:);DCN2(end-remove,:)]
C=[DCN3(remove+1,:);DCN3(end-remove,:)]

mns(1,:) = mean(DCN1(:,1:13));
mns(2,:) = mean(DCN2(:,1:13));
mns(3,:) = mean(DCN3(:,1:13));
lbs = [A(1,:)',B(1,:)',C(1,:)'];
ubs = [A(2,:)',B(2,:)',C(2,:)'];
s_table_diff_nc= strings([12,3]);
for i = 1:size(s_table_diff_nc,1)
    for j = 1:size(s_table_diff_nc,2)
        s_table_diff_nc(i,j) = strcat(num2str(mns(j,i),'%.3g'),{', ('},num2str(lbs(i,j),'%.3g'),{', '},num2str(ubs(i,j),'%.3g'),')');
    end
end

s_table_diff_nc
