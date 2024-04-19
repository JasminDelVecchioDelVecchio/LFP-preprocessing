function imu_embed = embedHeadAngularVelocity(vel, lags)
% input is a time x chan matrix

l = size(vel, 1); %time length
llag = length(lags); %

mi = min(lags);
ma = max(lags);

index = (-mi+1:l-ma)';
lindex = length(index);

imu_embed = [];
for ilag = 1:llag
    jj = index + lags(ilag);
  imu_embed = cat(2, imu_embed, vel(round(jj), :));
end

imu_embed = [imu_embed ones(lindex, 1)];
