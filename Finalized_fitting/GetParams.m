% generate 10,000 paramter estimates 

for k = 1:12
for j = 1:10000
idx = k
temp = ContactDataHaptoStimes(idx,1);
params(j,:,k) = temp(1,:);
params_noise(j,:,k) = temp(2,:);
REs(j,:,k) = temp(3,:);
temp = NeedleDataHaptoStimes(idx,1);
params_needle(j,:,k) = temp(1,:);
params_noise_needle(j,:,k) = temp(2,:);
REs_needle(j,:,k) = temp(3,:);
end
end

