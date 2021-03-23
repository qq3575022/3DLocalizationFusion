
clc, clear, close all

load('zyxCounClock/magx.mat'); load('zyxCounClock/magy.mat'); load('zyxCounClock/magz.mat');
load('zyxCounClock/gyrox.mat');load('zyxCounClock/gyroy.mat'); load('zyxCounClock/gyroz.mat');
load('zyxCounClock/accx.mat'); load('zyxCounClock/accy.mat');  load('zyxCounClock/accz.mat');
%%
accelReadings = zeros(length(accx), 3);
accelReadings(:, 1) = accx; accelReadings(:, 2) = accy; accelReadings(:, 3) = accz;


gyroReadings = zeros(length(gyrox), 3);
gyroReadings(:, 1) = gyrox;   gyroReadings(:, 2) = gyroy; gyroReadings(:, 3) = gyroz;


magReadings = zeros(length(magx), 3);
magReadings(:, 1) = magx; magReadings(:, 2) = magy; magReadings(:, 3) = magz;

FUSE = ahrsfilter

[orientation,angularVelocity] = FUSE(accelReadings,gyroReadings,magReadings);

figure
plot(orientation(1,:))