function [mag,gyro,AXY,Tmag,tdmag,Tgyro,tdgyro,Tacc,tdacc,T3,td3] = groundtruth3D2()

load('zyxClock/accT.mat');load('zyxClock/gyroT.mat');load('zyxClock/magT.mat');load('zyxClock/time.mat');

td3    = unique(round(time, 5));  %unique(round(time(time<2.1&time>0), 5)); 

tdacc  = round(accT,  5); %unique(round(accT(accT<2.1&accT>0), 5)); 
tdgyro = round(gyroT, 5); %unique(round(gyroT(gyroT<2.1&gyroT>0), 5)); 
tdmag  = round(magT,  5); %unique(round(magT(magT<2.1&magT>0), 5));

td3(1) = 0.00001; tdacc(1) = 0.00001; tdgyro(1) = 0.00001;
%

Tmag = zeros(1, length(tdmag)); Tacc = zeros(1, length(tdacc)); Tgyro = zeros(1, length(tdgyro)); T3 = zeros(1, length(td3)); 


for i = 2 : length(Tmag)
    Tmag(i)  = tdmag(i)  - tdmag(i-1);
end

for i = 2 : length(Tacc)
    Tacc(i)  = tdacc(i)  - tdacc(i-1);
end

for i = 2 : length(Tgyro)
    Tgyro(i) = tdgyro(i) - tdgyro(i-1);
end

for i = 2 : length(T3)
    T3(i)    = td3(i)    - td3(i-1);
end

load('zyxClock/magx.mat'); load('zyxClock/magy.mat'); load('zyxClock/magz.mat');
load('zyxClock/gyrox.mat');load('zyxClock/gyroy.mat'); load('zyxClock/gyroz.mat');
load('zyxClock/accx.mat'); load('zyxClock/accy.mat');  load('zyxClock/accz.mat');

mag = NaN(3,length(tdmag)); 
mag(1,:) = magx;   mag(2,:) = magy;   mag(3,:) = magz; 

gyro = NaN(3,length(tdgyro)); 
gyro(1,:) = gyrox; gyro(2,:) = gyroy;  gyro(3,:) = gyroz; 

AXY = NaN(3,length(tdacc)); 
AXY(1,:) = accx;   AXY(2,:) = accy;    AXY(3,:) = accz; 
end