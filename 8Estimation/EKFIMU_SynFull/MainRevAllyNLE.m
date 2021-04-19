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
yS = find(abs(acc_time-46.13)<0.002); yS = yS(1);yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);
yyS = find(abs(gyro_time-46.13)<0.002);  yyS = yyS(1);yyE = find(abs(gyro_time-48.5)<0.002); yyE = yyE(1);
ySS = find(abs(mag_time-46.13)<0.009);   ySS = ySS(1);yEE = find(abs(mag_time-48.5)<0.009); yEE = yEE(1);

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data)-g;
acc_data_z = data_z(find_acc_data);

gyro_data_x = data_x(find_gyro_data);gyro_data_y = data_y(find_gyro_data);gyro_data_z = data_z(find_gyro_data);
mag_data_x = data_x(find_mag_data);mag_data_y = data_y(find_mag_data);mag_data_z = data_z(find_mag_data);

% acc_data_x = medfilt1(acc_data_x,30);
% acc_data_y = medfilt1(acc_data_y,30);
% acc_data_z = medfilt1(acc_data_z,30);

accx = acc_data_x(yS:yE) - mean(acc_data_x(yS:yE));accy = acc_data_y(yS:yE) - mean(acc_data_y(yS:yE));accz = acc_data_z(yS:yE) - mean(acc_data_z(yS:yE));
accT = acc_time(yS:yE);

gyrox = gyro_data_x(yyS:yyE) - mean(gyro_data_x(yyS:yyE));gyroy = gyro_data_y(yyS:yyE) - mean(gyro_data_y(yyS:yyE));gyroz = gyro_data_z(yyS:yyE) - mean(gyro_data_z(yyS:yyE));
gyroT = gyro_time(yyS:yyE);


magx = mag_data_x(ySS:yEE) - mean(mag_data_x(ySS:yEE));magy = mag_data_y(ySS:yEE) - mean(mag_data_y(ySS:yEE));magz = mag_data_z(ySS:yEE) - mean(mag_data_z(ySS:yEE));
magT = mag_time(ySS:yEE);
% ========== Start End Time Acceleration ========== 

time = unique(sort([accT; gyroT; magT]),'rows');

[PP1, VV1, AA1] = groundtruth1Dy(time - time(1));

% Start from 0
PP1 = PP1 - PP1(1);
%
% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [xdotdot, ydotdot, orientation, angular velocity]
%x = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]'*ones(1,length(accT));
%o1 a1 o2 a2 o3 a3
x = zeros(15, length(time));
% Intial P matrix
P = 0.8*eye(length(x(:,1)));
% covariance for w
Q = diag([0.1, 0.1, 15, 0.1, 0.1, 15, 0.1, 0.1, 15, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]);

% covariance for v
RR = [0.3*var(accx), 0.3*var(accy), 0.3*var(accz), 0.03*var(magx), 0.3*var(gyrox), 0.03*var(magy), 0.3*var(gyroy), 0.03*var(magz), 0.3*var(gyroz)]; 

i = 1;j = 1;k = 1;
for m = 1:1:length(time)
    
    [y, i, j, k, index] = getyNPVA(magx, gyrox, magy, gyroy, magz, gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m);
    
    if m == 1
        x(:,m) = x(:,m);
    else
        F = getF(time, m);
 
        x(:,m) = F*x(:,m-1);

        P = F*P*F' + Q;
        H = getJacoN(y,x(:,m),index);

        R = getR(RR, index);

        G = P*H'/(H*P*H' + R);%lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m), index));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    
end

figure
subplot(311), plot(time, x(7,:), 'LineWidth', 3); hold on; plot(time, -PP1, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([46.13, 48.5]);
subplot(312), plot(time, x(8,:), 'LineWidth', 2); hold on; plot(time, -VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth','Location','SouthEast'); title('Velocity Along x Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([46.13, 48.5]);
subplot(313), plot(time, x(9,:), 'LineWidth', 2); hold on; plot(time, -AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth','Location','SouthEast'); title('Acceleration Along x Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([46.13, 48.5]);

rmsErrorY = rms(x(7,:)+PP1)
%
% figure
% subplot(9,1,1), plot(time, x(1,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
% subplot(9,1,2), plot(time, x(2,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
% subplot(9,1,3), plot(time, x(3,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(9,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
% subplot(9,1,4), plot(time, x(4,:), 'b', 'LineWidth', 2); %hold on; plot(time,, phi_gt(1,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
% subplot(9,1,5), plot(time, x(5,:), 'b', 'LineWidth', 2); %hold on; plot(time, phi_gt(2,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
% subplot(9,1,6), plot(time, x(6,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
% subplot(9,1,7), plot(time, x(7,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(4,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
% subplot(9,1,8), plot(time, x(8,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(5,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
% subplot(9,1,9), plot(time, x(9,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;
