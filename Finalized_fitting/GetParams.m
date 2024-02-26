% generate 10,000 paramter estimates 
% loop over needle and contact infected hosts 
for k = 1:12
for j = 1:10000
idx = k
% estimate for each contact invected host 
temp = ContactDataHaptoStimes(idx,1);
% store baseline estimate 
params(j,:,k) = temp(1,:);
% store estimate from generated dataset reestimation for practicial
% identifiability and interval estimates 

% do the same for each needle infected host 
params_noise(j,:,k) = temp(2,:);
REs(j,:,k) = temp(3,:);
temp = NeedleDataHaptoStimes(idx,1);
params_needle(j,:,k) = temp(1,:);
params_noise_needle(j,:,k) = temp(2,:);
end
end

