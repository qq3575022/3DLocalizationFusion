clc, clear, close all
data=readtable('3.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));               
data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

% start time
xS = find(abs(acc_time-3.5)<0.002); xS = xS(1); xE = find(abs(acc_time-7.5)<0.001); xE = xE(1);
xxS = find(abs(gyro_time-3.5)<0.002); xxS = xxS(1); xxE = find(abs(gyro_time-7.5)<0.001); xxE = xxE(1);
xSS = find(abs(mag_time-3.5)<0.002); xSS = xSS(1); xEE = find(abs(mag_time-7.5)<0.004); xEE = xEE(1);

yS = find(abs(acc_time-46)<0.002); yS = yS(1); yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);
zS = find(abs(acc_time-95.5)<0.002); zS = zS(1); zE = find(abs(acc_time-98)<0.002); zE = zE(1);

%time = unique(sort([accT; gyroT; magT]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data)-g;
acc_data_z = data_z(find_acc_data);

gyro_data_x = data_x(find_gyro_data); gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

accx = acc_data_x(xS:xE) - mean(acc_data_x(xS:xE));accy = acc_data_y(xS:xE) - mean(acc_data_y(xS:xE));accz = acc_data_z(xS:xE) - mean(acc_data_z(xS:xE));
accT = acc_time(xS:xE);

gyrox = gyro_data_x(xxS:xxE) - mean(gyro_data_x(xxS:xxE));gyroy = gyro_data_y(xxS:xxE) - mean(gyro_data_y(xxS:xxE));gyroz = gyro_data_z(xxS:xxE) - mean(gyro_data_z(xxS:xxE));
gyroT = gyro_time(xxS:xxE);

magx = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE));magy = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE));magz = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));
magT = mag_time(xSS:xEE);
% ========== Start End Time Acceleration ========== 

time = unique(sort([accT; gyroT; magT]),'rows');
%
[PP1, VV1, AA1] = groundtruth1Dx(time - time(1)); PP1 = PP1 - PP1(1);

N = 5; i = 1; j = 1; k = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);
    
for m = 1:1:length(time)-N
    
    [y, i, j, k] = gety(magx, gyrox, magy, gyroy, magz, gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, time, m, accT, i),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+1);
    x(:,m) = F*x(:,m);

end
%%
[pos, vel] = doubleInte(diff(time(1:end-N)), x(3,2:end));


figure
subplot(311), plot(time(2:end-N), pos, 'LineWidth', 3); hold on; plot(time, PP1, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','SouthEast'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([3.5, 7.5])
subplot(312), plot(time(2:end-N), vel, 'LineWidth', 2); hold on; plot(time, VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth'); title('Velocity Along x Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([3.5, 7.5])
subplot(313), plot(time(1:end-N), x(3,1:end), 'LineWidth', 2); hold on; plot(time, AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth'); title('Acceleration Along x Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([3.5, 7.5])

rmsErrorX = rms(pos-PP1(2:end-N))
%
figure
subplot(9,1,1), plot(time(1:end-N), x(1,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
subplot(9,1,2), plot(time(1:end-N), x(2,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
subplot(9,1,3), plot(time(1:end-N), x(3,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(9,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
subplot(9,1,4), plot(time(1:end-N), x(4,:), 'b', 'LineWidth', 2); %hold on; plot(time,, phi_gt(1,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
subplot(9,1,5), plot(time(1:end-N), x(5,:), 'b', 'LineWidth', 2); %hold on; plot(time, phi_gt(2,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
subplot(9,1,6), plot(time(1:end-N), x(6,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
subplot(9,1,7), plot(time(1:end-N), x(7,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(4,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
subplot(9,1,8), plot(time(1:end-N), x(8,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(5,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
subplot(9,1,9), plot(time(1:end-N), x(9,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;
