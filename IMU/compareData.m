clc, clear, close all
% import csv data and type in gravity
data=readtable('../Data/0304IMU/whole/2first.csv','Delimiter', ',');

g=9.7953;%https://www.sensorsone.com/local-gravity-calculator/

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)_(Axis Pointing Up When Test Started/ Axis being investigated)
%Extract data when z axis bing investigated, name by acc_data_(Axis)_zt
sensor_acc_check = ismember(data.Var2,'ACC_UN');
find_acc_data_zt = find(sensor_acc_check == 1);

data_time = table2array(data(:,1));
data_x_zt = table2array(data(:,3));
data_y_zt = table2array(data(:,4));
data_z_zt = table2array(data(:,5));

data_time = data_time(find_acc_data_zt);
accT = (data_time - data_time(1))/1000000000;

acc_data_x_zt=data_x_zt(find_acc_data_zt);
acc_data_y_zt=data_y_zt(find_acc_data_zt);
acc_data_z_zt=data_z_zt(find_acc_data_zt);


% figure;
% plot(accT, acc_data_x_zt)
% hold on
% plot(accT, acc_data_y_zt)
% hold on
% plot(accT, acc_data_z_zt)
% xlabel('samples')
% ylabel('sensor read')
% legend('x-aixs','y-aixs','z-aixs')
% title('accelerometer data, z-axis fliped')

accx = medfilt1(acc_data_x_zt,40);
accy = medfilt1(acc_data_y_zt,40);
accz = medfilt1(acc_data_z_zt,40);

% accx = accx - accx(10) + 0.005 - 0.000003*accT;
% accy = accy - accy(10) - 0.018 - 0.00004*accT;
% accz = accz - accz(10) - 0.065 + 0.00005*accT;

driftX = accx(find(abs(accT-35)<0.001)) - accx(5)
driftX/35

driftY = accy(find(abs(accT-35)<0.001)) - accy(5)
driftY/35

driftZ = accz(find(abs(accT-35)<0.001)) - accz(5)
driftZ/35

% posX = doubleInte(diff(accT), accx(2:length(accx)), 2);
% posY = doubleInte(diff(accT), accy(2:length(accy)), 2);
% posZ = doubleInte(diff(accT), accz(2:length(accz)), 2);
% 
% figure
% subplot(311), plot(accT(2:length(accT)), posX)
% subplot(312), plot(accT(2:length(accT)), posZ)
% subplot(313), plot(accT(2:length(accT)), posY)
% figure
% plot(posX, posZ)

figure;
plot(accT, acc_data_x_zt)
hold on
plot(accT, acc_data_y_zt)
hold on
plot(accT, acc_data_z_zt)
hold on
plot(accT, accx)
hold on
plot(accT, accy)
hold on
plot(accT, accz)
xlabel('samples')
ylabel('sensor read')
legend('x-aixs','y-aixs','z-aixs','x-aixs','y-aixs','z-aixs')
title('accelerometer data, z-axis fliped')
% 
% save('accx.mat','accx');
% save('accy.mat','accy');
% save('accz.mat','accz');
% 
% save('accT.mat','accT');


%%
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');
find_gyro_data_zt = find(sensor_gyro_check == 1);

data_time = table2array(data(:,1));
data_x_gyro = table2array(data(:,3));
data_y_gyro = table2array(data(:,4));
data_z_gyro = table2array(data(:,5));

data_time = data_time(find_gyro_data_zt);
gyroT = (data_time - data_time(1))/1000000000;

gyro_data_x=data_x_gyro(find_gyro_data_zt);
gyro_data_y=data_y_gyro(find_gyro_data_zt);
gyro_data_z=data_z_gyro(find_gyro_data_zt);


% figure;
% plot(gyroT, gyro_data_x)
% hold on
% plot(gyroT, gyro_data_y)
% hold on
% plot(gyroT, gyro_data_z)
% xlabel('samples')
% ylabel('sensor read')
% legend('x-aixs','y-aixs','z-aixs')
% title('accelerometer data, z-axis fliped')

gyrox = medfilt1(gyro_data_x,40);
gyroy = medfilt1(gyro_data_y,40);
gyroz = medfilt1(gyro_data_z,120);

% x = averaging_filter_mex(acc_data_x_zt);
% y = averaging_filter_mex(acc_data_y_zt);
% z = averaging_filter_mex(acc_data_z_zt);

% [timex, x3, x4, x5, x6, x7] = filterAvg(data_x_zt);
% [timey, y3, y4, y5, y6, y7] = filterAvg(data_y_zt);
% [timez, z3, z4, z5, z6, z7] = filterAvg(data_z_zt);


% figure;
% subplot(311), plot(gyroT, gyro_data_x), ylim([-0.15, 0.15]);
% hold on
% plot(gyroT, gyrox)
% xlabel('samples')
% ylabel('sensor read')
% legend('x-aixs','filtered')
% title('accelerometer data, x-axis')
% 
% subplot(312), plot(gyroT, gyro_data_y), ylim([-0.15, 0.15]);
% hold on
% plot(gyroT, gyroy)
% xlabel('samples')
% ylabel('sensor read')
% legend('y-aixs','filtered')
% title('accelerometer data, y-axis')
% 
% subplot(313), plot(gyroT, gyro_data_z), ylim([-0.15, 0.15]);
% hold on
% plot(gyroT, gyroz)
% xlabel('samples')
% ylabel('sensor read')
% legend('z-aixs','filtered')
% title('accelerometer data, z-axis')
% 
% save('gyrox.mat','gyrox');
% save('gyroy.mat','gyroy');
% save('gyroz.mat','gyroz');
% 
% save('gyroT.mat','gyroT');

%load('xyzClock/accT.mat');load('xyzClock/gyroT.mat');load('xyzClock/magT.mat');load('xyzClock/time.mat');


%% Mag

sensor_gyro_check = ismember(data.Var2,'MAG_UN');
find_gyro_data_zt = find(sensor_gyro_check == 1);
data_time = table2array(data(:,1));
data_x_gyro = table2array(data(:,3));
data_y_gyro = table2array(data(:,4));
data_z_gyro = table2array(data(:,5));

data_time = data_time(find_gyro_data_zt);
data_time = (data_time - data_time(1))/1000000000;

gyro_data_x=data_x_gyro(find_gyro_data_zt);
gyro_data_y=data_y_gyro(find_gyro_data_zt);
gyro_data_z=data_z_gyro(find_gyro_data_zt);


figure;
plot(data_time, gyro_data_x)
hold on
plot(data_time, gyro_data_y)
hold on
plot(data_time, gyro_data_z)
xlabel('samples')
ylabel('sensor read')
legend('x-aixs','y-aixs','z-aixs')
title('accelerometer data, z-axis fliped')

x = medfilt1(gyro_data_x,40);
y = medfilt1(gyro_data_y,40);
z = medfilt1(gyro_data_z,120);

% x = averaging_filter_mex(acc_data_x_zt);
% y = averaging_filter_mex(acc_data_y_zt);
% z = averaging_filter_mex(acc_data_z_zt);

% [timex, x3, x4, x5, x6, x7] = filterAvg(data_x_zt);
% 
% [timey, y3, y4, y5, y6, y7] = filterAvg(data_y_zt);
% 
% [timez, z3, z4, z5, z6, z7] = filterAvg(data_z_zt);


% figure;
% subplot(311), plot(data_time, gyro_data_x), ylim([-0.15, 0.15]);
% hold on
% plot(data_time, x)
% 
% subplot(312), plot(data_time, gyro_data_y), ylim([-0.15, 0.15]);
% hold on
% plot(data_time, y)
% 
% subplot(313), plot(data_time, gyro_data_z), ylim([-0.15, 0.15]);
% hold on
% plot(data_time, z)
% xlabel('samples')
% ylabel('sensor read')
% legend('x-aixs','y-aixs','z-aixs','x-aixs','y-aixs','z-aixs')
% title('accelerometer data, z-axis fliped')

