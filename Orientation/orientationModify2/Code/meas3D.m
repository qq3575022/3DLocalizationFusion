function [gyro,AXY,Tgyro,tdgyro,Tacc,tdacc] = groundtruth3D2()


%data=readtable('../3DLocalizationFusion/Data/0304IMU/whole/1.csv','Delimiter', ',');  g=9.7953;
data=readtable('1.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

tdacc   = table2array(data(:,1));         data_x     = table2array(data(:,3)); 
tdacc   = tdacc(find_acc_data);           data_y     = table2array(data(:,4));    
tdacc   = (tdacc - tdacc(1))/1000000000;  data_z     = table2array(data(:,5));

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
accy = medfilt1(acc_data_y,500);
accz = medfilt1(acc_data_z,500);
%%
% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  

tdgyro   = table2array(data(:,1));         data_x     = table2array(data(:,3)); 
tdgyro   = tdgyro(find_acc_data);           data_y     = table2array(data(:,4));    
tdgyro   = (tdgyro - tdgyro(1))/1000000000;  data_z     = table2array(data(:,5));

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
gyro_data_x = data_x(find_gyro_data);
gyro_data_y = data_y(find_gyro_data);
gyro_data_z = data_z(find_gyro_data);

% ========== Filtered Acceleration Data ========== 
% ----  accx  ----
% ----  accy  ----
% ----  accz  ----
gyrox = medfilt1(gyro_data_x,500);
gyroy = medfilt1(gyro_data_y,500);
gyroz = medfilt1(gyro_data_z,500);
%%
tdacc  = round(tdacc,  5); %unique(round(accT(accT<2.1&accT>0), 5)); 
tdgyro = round(tdgyro, 5); %unique(round(gyroT(gyroT<2.1&gyroT>0), 5)); 
%tdmag  = round(magT,  5); %unique(round(magT(magT<2.1&magT>0), 5));

%

%Tmag = zeros(1, length(tdmag)); T3 = zeros(1, length(td3)); 
Tacc = zeros(1, length(tdacc)); Tgyro = zeros(1, length(tdgyro)); 


% for i = 2 : length(Tmag)
%     Tmag(i)  = tdmag(i)  - tdmag(i-1);
% end

for i = 2 : length(Tacc)
    Tacc(i)  = tdacc(i)  - tdacc(i-1);
end

for i = 2 : length(Tgyro)
    Tgyro(i) = tdgyro(i) - tdgyro(i-1);
end



% mag = NaN(3,length(tdmag)); 
% mag(1,:) = magx;   mag(2,:) = magy;   mag(3,:) = magz; 

gyro = NaN(3,length(tdgyro)); 
gyro(1,:) = gyrox; gyro(2,:) = gyroy;  gyro(3,:) = gyroz; 

AXY = NaN(3,length(tdacc)); 
AXY(1,:) = accx;   AXY(2,:) = accy;    AXY(3,:) = accz; 
end