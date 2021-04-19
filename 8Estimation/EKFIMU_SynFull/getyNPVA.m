function [y, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m)
% m = m
% i = i
% j = j
% k = k
% timeTt = time(m)
% accTt = acc_time(i)
% gyroTt = gyro_time(j)
% magTt = mag_time(k)

if i == 1
    indexI = 1;
elseif i <= length(acc_time) && acc_time(i-1) >= time(m)
    indexI = i-1;
elseif i <= length(acc_time) && acc_time(i-1) < time(m)
    indexI = i;
else
    indexI = i-1;
end

if j == 1
    indexJ = 1;
elseif j <= length(gyro_time) && gyro_time(j-1) >= time(m)
    indexJ = j-1;
elseif j <= length(gyro_time) && gyro_time(j-1) < time(m)
    indexJ = j;
else
    indexJ = j-1;
end

if k == 1
    indexK = 1;
elseif k <= length(mag_time) && mag_time(k-1) >= time(m)
    indexK = k-1;
elseif k <= length(mag_time) && mag_time(k-1) < time(m)
    indexK = k;
else
    indexK = k-1;
end


%if i <= length(acc_time) && j <= length(gyro_time) && k <= length(mag_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - gyro_time(j)) ==0 && abs(time(m) - mag_time(k)) ==0
y = NaN(9,1);

y(1:9) = [mag_data_x(indexK), gyro_data_x(indexJ), mag_data_y(indexK), gyro_data_y(indexJ), mag_data_z(indexK), gyro_data_z(indexJ), acc_data_x(indexI), acc_data_y(indexI), acc_data_z(indexI)];

i = indexI + 1;
j = indexJ + 1;
k = indexK + 1;

end