clc, clear, close all
data=readtable('../Data/0304IMU/whole/1.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

acc_time   = table2array(data(:,1));               data_x     = table2array(data(:,3)); 
acc_time   = acc_time(find_acc_data);              data_y     = table2array(data(:,4));    
acc_time   = (acc_time - acc_time(1))/1000000000;  data_z     = table2array(data(:,5));

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data);
acc_data_z = data_z(find_acc_data);

% ========== Filtered Acceleration Data ========== 
% ----  accx  ----
% ----  accy  ----
% ----  accz  ----
accx = medfilt1(acc_data_x,500);
accy = medfilt1(acc_data_y,500) - g;
accz = medfilt1(acc_data_z,500);

%
biasx = 0.014;     biasy = 0.028;    biasz = 0.075;
scalarx = 1;   scalary = 1;  scalarz = 1; 
skew12  = 0;   skew13  = 0;  skew23  = 0; 
%driftx = 0.002; drifty = 0.001; driftz = 0.001;
driftx = 0.00; drifty = 0.00; driftz = 0.00;

cali = [scalarx, skew12, skew13; 0, scalary, skew23; 0, 0, scalarz]*[(acc_data_x-mean(acc_data_x)-biasx)'; (acc_data_y-mean(acc_data_y)-biasy)'; (acc_data_z-mean(acc_data_z)-biasz)'];

accxCali = cali(1, :) + driftx*acc_time';
accyCali = cali(2, :) + drifty*acc_time';
acczCali = cali(3, :) + driftz*acc_time';
% ******************************  //1 Comment out // } ******************************


driftxSep = 0.078; driftySep = 0.03; driftzSep = 0.055; 


% ========== Start End Time Acceleration ========== 
xS = find(abs(acc_time-35)<0.002); xS = xS(1);
xE = find(abs(acc_time-48)<0.001); xE = xE(1);

yS = find(abs(acc_time-50.5178)<0.002); yS = yS(1);
yE = find(abs(acc_time-68.5178)<0.002); yE = yE(1);

zS = find(abs(acc_time-70.5178)<0.002); zS = zS(1);
zE = find(abs(acc_time-83.5178)<0.002); zE = zE(1);

BS = find(abs(acc_time-85.5178)<0.002); BS = BS(1);
BE = find(abs(acc_time-98.5178)<0.002); BE = BE(1);
BEY =find(abs(acc_time-105.5178)<0.002);BEY = BEY(1);

% ========== Calibrated Acceleration Data ========== 
% ---- accxCali  ----
% ---- accyCali  ----
% ---- acczCali  ----

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Method 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% unknow = [0.98, 0.02, 0.01; 0, 0.92, 0.02; 0, 0, 0.96]*[(accx - 0.02)'; (accy + 0.01)'; (accz  - 0.03)'];
% ---  biasx biasy biasz --- 
% --- driftx drifty driftz ---

% driftX = accx(find(abs(accT-35)<0.001)) - accx(5)
% driftX/35
% 
% driftY = accy(find(abs(accT-35)<0.001)) - accy(5)
% driftY/35
% 
% driftZ = accz(find(abs(accT-35)<0.001)) - accz(5)
% driftZ/35

% ******************************{ //1 Comment out // ******************************
biasx = 0.014;     biasy = 0.028;    biasz = 0.075;
scalarx = 1;   scalary = 1;  scalarz = 1; 
skew12  = 0;   skew13  = 0;  skew23  = 0; 
%driftx = 0.002; drifty = 0.001; driftz = 0.001;
driftx = 0.00; drifty = 0.00; driftz = 0.00;

cali = [scalarx, skew12, skew13; 0, scalary, skew23; 0, 0, scalarz]*[(accx-mean(accx)-biasx)'; (accy-mean(accy)-biasy)'; (accz-mean(accz)-biasz)'];

accxCali = cali(1, :) + driftx*acc_time';
accyCali = cali(2, :) + drifty*acc_time';
acczCali = cali(3, :) + driftz*acc_time';
% ******************************  //1 Comment out // } ******************************


driftxSep = 0.078; driftySep = 0.03; driftzSep = 0.055; 

% ============= Integration to Position ============= 
timeDiff = diff(acc_time);
posX  = doubleInte(timeDiff(xS:xE), accxCali(xS:xE), 2);
posX2 = doubleInte(timeDiff(BS:BE), accxCali(BS:BE) + driftxSep, 2);

posY  = doubleInte(timeDiff(yS:yE), acczCali(yS:yE), 2);
posY2 = doubleInte(timeDiff(BS:BE), acczCali(BS:BE) + driftzSep, 2);

posZ  = doubleInte(timeDiff(zS:zE), accyCali(zS:zE), 2);
posZ2 = doubleInte(timeDiff(BS:BE), accyCali(BS:BE) + driftySep, 2);


posYFinal  = doubleInte(timeDiff(yS:BEY), acczCali(yS:BEY)-0.015, 2);

figure
plot(posYFinal)

figure
subplot(611), plot(acc_time(xS:xE), posX, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('x move'), grid on, grid minor
subplot(612), plot(acc_time(BS:BE), posX2, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('x back'), grid on, grid minor

subplot(613), plot(acc_time(yS:yE), posY, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('y move'), grid on, grid minor
subplot(614), plot(acc_time(BS:BE), posY2, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('y back'), grid on, grid minor
 

subplot(615), plot(acc_time(zS:zE), posZ, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('z move'), grid on, grid minor
subplot(616), plot(acc_time(BS:BE), posZ2, 'LineWidth', 2), ylabel('position [m]'), xlabel('time [s]'), title('z back'), grid on, grid minor


figure;
plot(acc_time(2:length(acc_time)), accx(2:length(accx)), 'LineWidth', 2)
hold on
plot(acc_time(2:length(acc_time)), accy(2:length(accx)), 'LineWidth', 2)
hold on
plot(acc_time(2:length(acc_time)), accz(2:length(accx)), 'LineWidth', 2)
xlabel('samples')
ylabel('sensor read')
legend('x-aixs','y-aixs','z-aixs')
title('accelerometer data')

% save('accx.mat','accx');
% save('accy.mat','accy');
% save('accz.mat','accz');
% 
% save('accT.mat','accT');

% figure
% plot(posX, posZ)
% figure;
% plot(accT, acc_data_x_zt)
% hold on
% plot(accT, acc_data_y_zt)
% hold on
% plot(accT, acc_data_z_zt)
% hold on
% plot(accT, accx)
% hold on
% plot(accT, accy)
% hold on
% plot(accT, accz)
% xlabel('samples')
% ylabel('sensor read')
% legend('x-aixs','y-aixs','z-aixs','x-aixs','y-aixs','z-aixs')
% title('accelerometer data, z-axis fliped')
