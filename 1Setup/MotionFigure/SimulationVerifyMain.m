clc, clear, close all

load('PropGp1DVelocityAccelerationData.mat')
ta = 0.4657;
tb = 1.1657;

td = Time(Time >= ta & Time <= tb)'-ta;
T = mean(diff(td));

td  = 0:T:13;
td1 = 0:T:23;
td3 = 0:T:13*2+23*2;

[PP1, VV1, AA1] = groundtruth1DxRF(td);
[PP2, VV2, AA2] = groundtruth1DyRF(td1);
[PP3, VV3, AA3] = groundtruth1DzRF(td);

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

%
% 1D - 3D coordinates
% td3
A3 =  zeros(3, 2*length(td)+2*length(td1));

% x
A3(1,1:length(AA1)) = AA1;           %      A3(1,length(AA1)*2+length(AA2)+1:length(AA1)*3+length(AA2)) = flip(AA1); 
%y
A3(2,length(AA1)+1:length(AA1) + length(AA2)) = AA2;  % A3(2,length(AA1)*2 + length(AA2)+1:end) = flip(AA2);
%z
A3(3,length(AA1) + length(AA2)+1:length(AA1)*2 + length(AA2)) = AA3;% A3(2,length(AA1)*2+length(AA2)+1:length(AA1)*3+length(AA2)) = flip(AA3);

%

% 3D coordinates
coord3 = zeros(3, length(td3));
%
% x
coord3(1,1:length(AA1)) = PP1 + 1.03;
coord3(1,length(AA1)+1:length(AA1)+length(AA2)) = PP1(end) + 1.03;
coord3(1,length(AA1)+length(AA2)+1:length(AA1)*2+length(AA2)) = PP1(end) + 1.03;%flip(PP1) + 1.03;
coord3(1,length(AA1)*2+length(AA2)+1:end) = PP1(end) + 1.03;%1.03;
%y
coord3(2,1:length(AA1)) = 1.31;
coord3(2,length(AA1)+1:length(AA1)+length(AA2)) = PP2 + 1.31;
coord3(2,length(AA1)+length(AA2)+1:length(AA1)*2+length(AA2)) = PP2(end) + 1.31;
coord3(2,length(AA1)*2+length(AA2)+1:end) = PP2(end) + 1.31;%flip(PP2) + 1.31;
%z
coord3(3,1:length(AA1)+length(AA2)) = 1.03;
coord3(3,length(AA1)+length(AA2)+1:length(AA1)*2+length(AA2)) = PP3 + 1.03;
coord3(3,length(AA1)*2+length(AA2)+1:end) = PP3(end) + 1.03;%flip(PP3) + 1.03;
%coord3(3,length(AA1)*3+length(AA2)+1:end) = PP3(end) + 1.03;%1.03;


% Dash Dot
coord3_ = zeros(3, length(AA2));
coord3_(1,1:length(AA1)) = flip(PP1) + 1.03;
coord3_(1,length(AA1)+1:end) = 1.03;
%y
coord3_(2,1:length(AA2)) = flip(PP2) + 1.31;
%z
coord3_(3,1:length(AA1)) = flip(PP3) + 1.03;
coord3_(3,length(AA1)+1:end) = 1.03;

%
% position of readers
%     x1 = [-0.05, 1.5];x2 = [2, 3];x3 = [2.7, 0.05];
x1 = [0,    0,     0.865];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,     1.27];%[ 0.9, 0.8, 0.2];
x3 = [2.29, 2.52,  0.865];%[ 0.8,-0.3,-0.2];
x4 = [0,    2.52,  1.27];
% x4 = [-0.3, 0.9, 0.8];


figure
plot3(coord3(1,:), coord3(2,:), coord3(3,:), 'LineWidth', 3);xlim([0 2.5]);ylim([0 2.7]);zlim([0.85 1.28])
hold on;
plot3(coord3_(1,:), coord3_(2,:), coord3_(3,:), 'LineWidth', 1, 'color', [0, 0.4470, 0.7410], "LineStyle",'-.');
hold on;
scatter3(x1(1), x1(2), x1(3), '*'); text(x1(1), x1(2) + 0.1, x1(3),"Reader 1 (0,0, 0.865)",'fontsize',16);
hold on;
scatter3(x2(1), x2(2), x2(3), '*'); text(x2(1), x2(2) + 0.1, x2(3),"Reader 2 (2.29,0,1.27)",'fontsize',16);
hold on;
scatter3(x3(1), x3(2), x3(3), '*'); text(x3(1) - 0.7, x3(2) - 0.5, x3(3) + 0.03,"Reader 3 (2.29,2.52,0.865)",'fontsize',16);
hold on;
scatter3(x4(1), x4(2), x4(3), '*'); text(x4(1), x4(2) + 0.1, x4(3),"Reader 4 (0,2.52,1.27)",'fontsize',16);
grid on;
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
%
%save('coord3.mat','coord3')
% figure
% subplot(311), plot(td3, coord3(1,:))
% subplot(312), plot(td3, coord3(2,:))
% subplot(313), plot(td3, coord3(3,:))