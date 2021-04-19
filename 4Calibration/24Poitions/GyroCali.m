clc, clear, close all

load('24Winpos4.mat')
load('24Gyrpos4.mat')
load('24gTpos4.mat')
load('24Stpos4.mat')
% 
data=readtable('24pos4.csv','Delimiter', ',');

sensor_gyro_check = ismember(data.Var2,'GYRO_UN');
find_gyro_data_zt = find(sensor_gyro_check == 1);

data_z = table2array(data(:,5));
data_y = table2array(data(:,4));
data_x = table2array(data(:,3));

gyro_z=data_z(find_gyro_data_zt);
gyro_y=data_y(find_gyro_data_zt);
gyro_x=data_x(find_gyro_data_zt);

% cali res

A = [0.9875, 0.0825, 0.0497, 0.1879;0, 0.9964, -0.0481, 0.0555; 0, 0, 0.9845, -0.2620];
y = [res ones(length(res),1)];

yy = A * y';

x=lsqnonlin(@(xx)getEG(yy, xx, res, resgT, gyro_x, gyro_y, gyro_z, windows),[0.1;0.1;0.1;0.1;0.1;0.1;1;1;1]);

%%
gCali = [x(7), x(1)*x(8), x(2)*x(9); x(3)*x(7), x(8), x(4)*x(9); x(5)*x(7), x(6)*x(8), x(9)];