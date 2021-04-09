clc, clear, close all
% get angular velocity w
load('X.mat');load('A.mat');load('AxB.mat');load('AyB.mat')
load('yVectorData.mat')
%%
h = 1e-5;           tc = 0:h:2.00000;
T = mean(diff(tf)); td = 0:T:2.00000;

w_a = zeros(1,length(tc));w_v = NaN(1,length(tc));w = NaN(1,length(tc));
W_A = NaN(1,length(td));  W_V = NaN(1,length(td));W = NaN(1,length(td));
w(1) = 0; w_v(1) = 0;

w_a(tc > 0.16 & tc <= 0.74) = 12.5;
w_a(tc > 0.74 & tc <= 1.05) = 0;
w_a(tc > 1.05 & tc <= 1.63) = -12.5;

figure
plot(tc,A)
hold on;
plot(tc,w_a)

% w_a(tc > 0.13 & tc <= 0.74) = 10.9;
% w_a(tc > 0.74 & tc <= 1.08) = 0;
% w_a(tc > 1.08 & tc <= 1.69) = -10.9;


