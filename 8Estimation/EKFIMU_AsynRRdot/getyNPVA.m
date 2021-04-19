function [y, i, j, k, index, len] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m)
% m = m
% i = i
% j = j
% k = k
% timeTt = time(m)
% accTt = acc_time(i)
% gyroTt = gyro_time(j)
% magTt = mag_time(k)

if i <= length(acc_time) && j <= length(gyro_time) && k <= length(mag_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - gyro_time(j)) ==0 && abs(time(m) - mag_time(k)) ==0
  index = 1;
  len = 9;
  y = NaN(9,1);
  y(1:9) = [mag_data_x(k), gyro_data_x(j), mag_data_y(k), gyro_data_y(j), mag_data_z(k), gyro_data_z(j), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
  i = i + 1;
  j = j + 1;
  k = k + 1;
  
elseif i <= length(acc_time) && j <= length(gyro_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - gyro_time(j)) == 0
  index = 2;
  len = 6;
  y = NaN(6,1);
  y(1:6) = [gyro_data_x(j), gyro_data_y(j), gyro_data_z(j), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
  i = i + 1;
  j = j + 1;
  
elseif j <= length(gyro_time) && k <= length(mag_time) && abs(time(m) - gyro_time(j)) == 0 && abs(time(m) - mag_time(k)) == 0
  index = 3;
  len = 6;
  y = NaN(6,1);
  y(1:6) = [mag_data_x(k), gyro_data_x(j), mag_data_y(k), gyro_data_y(j), mag_data_z(k), gyro_data_z(j)];
  j = j + 1;
  k = k + 1;
  
elseif i <= length(acc_time) && k <= length(mag_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - mag_time(k)) == 0
  index = 4;
  len = 6;
  y = NaN(6,1);
  y(1:6) = [mag_data_x(k), mag_data_y(k), mag_data_z(k), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
  i = i + 1;
  k = k + 1;
  
elseif i <= length(acc_time) && abs(time(m) - acc_time(i)) == 0
  index = 5;
  len = 3;
  y = NaN(3,1);
  y(1:3) = [acc_data_x(i), acc_data_y(i), acc_data_z(i)];
  i = i + 1;
  
elseif j <= length(gyro_time) && abs(time(m) - gyro_time(j)) < 0.0001
  index = 6;
  len = 3;
  y = NaN(3,1);
  y(1:3) = [gyro_data_x(j), gyro_data_y(j), gyro_data_z(j)];
  j = j + 1;
  
elseif k <= length(mag_time) && abs(time(m) - mag_time(k)) < 0.0001
  index = 7;
  len = 3;
  y = NaN(3,1);
  y(1:3) = [mag_data_x(k), mag_data_y(k), mag_data_z(k)];
  k = k + 1;
  
else
    y = nan;
    len = -1;
    index = -1;
  
end

end