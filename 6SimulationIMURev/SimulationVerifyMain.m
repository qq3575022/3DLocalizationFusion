clc, clear, close all

% compare
data=readtable('3.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

acc_time   = table2array(data(:,1));               data_x     = table2array(data(:,3)); 
acc_time   = acc_time(find_acc_data);              data_y     = table2array(data(:,4));    
acc_time   = (acc_time - acc_time(1))/1000000000;  data_z     = table2array(data(:,5));

% 3.5  7.5   46  48.5   95.5  98
xS = find(abs(acc_time-3.6)<0.002); xS = xS(1);
xE = find(abs(acc_time-7.5)<0.001); xE = xE(1);

yS = find(abs(acc_time-46.1784)<0.002); yS = yS(1);
yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);

zS = find(abs(acc_time-95.5027)<0.002); zS = zS(1);
zE = find(abs(acc_time-98)<0.002); zE = zE(1);


% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data);
acc_data_z = data_z(find_acc_data);

%
%td = 0:T:5;
tdx = acc_time(xS:xE) - acc_time(xS);
tdy = acc_time(yS:yE) - acc_time(yS);
tdz = acc_time(zS:zE) - acc_time(zS);
%%
[PP1, VV1, AA1] = groundtruth1Dx(tdx);
[PP2, VV2, AA2] = groundtruth1Dy(tdy);
[PP3, VV3, AA3] = groundtruth1Dz(tdz);

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

% 1D - 3D coordinates
% td3

A3 =  zeros(3, length(acc_time));

%xS = 3.6; yS = 46.1784; zS = 95.5027;

% x
A3(1,xS+1:length(AA1) + xS) = AA1;
%y
A3(2,yS+1:length(AA2) + yS) = AA2;
%z
A3(3,zS+1:length(AA3) + zS) = AA3;
% A3(3,length(AA3)*3+1:length(AA3)*4) = flip(AA3);

%     Ta = [1 0.0936 0.0621;
%           0 1      0.0329;
%           0 0      1];
%     Ka =  diag([1.001, 0.9812, 0.9441]);
%     
%     TKa = pinv(Ta*Ka);
%     
%     errorA = TKa*A3;
%     errorA3 = errorA - [0.0961;0.0928;0.0947];



% ========== Filtered Acceleration Data ========== 
% ----  accx  ----
% ----  accy  ----
% ----  accz  ----
accx = medfilt1(acc_data_x,50);
accy = medfilt1(acc_data_y,50) - g;
accz = medfilt1(acc_data_z,50);

% xS = find(abs(acc_time-5)<0.002); xS = xS(1);
% BEY =find(abs(acc_time-90.5178)<0.002);BEY = BEY(1);

A3 = A3 + [0.14273;0.650848;0]*ones(1, length(acc_time));
A3_sim = A3 + [0;-0.0005;0]*acc_time';

A3_sim = awgn(A3_sim,35,'measured');

figure
subplot(611),plot(acc_time, A3_sim(1,:),'LineWidth',2), xlim([0,100]);%hold on, plot(acc_time, accx,'LineWidth',2), xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(612),plot(acc_time, accx,'LineWidth',2),xlim([0,100]);%hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(613),plot(acc_time, A3_sim(2,:),'LineWidth',2), %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(614),plot(acc_time, accz,'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(615),plot(acc_time, A3_sim(3,:),'LineWidth',2), %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')
subplot(616),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')

figure
subplot(311),plot(acc_time, accx,'LineWidth',3),hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(312),plot(acc_time, accz,'LineWidth',3),hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(313),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',3),hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')

%%

% 1D - 3D coordinates
coord3 = zeros(3, length(acc_time));

% x
coord3(1,round(3.28745/T)+1:length(AA1) + round(3.28745/T)) = PP1 + 1.03;
coord3(1,length(AA1) + round(3.28745/T) + 1:end) = PP1(end) + 1.03;

%y
coord3(2,1:round(45.8019/T)) = 1.31;
coord3(2,round(45.8019/T)+1:length(AA2) + round(45.8019/T)) = PP2 + 1.31;
coord3(2,length(AA2) + round(45.8019/T)+1:end) = PP2(end) + 1.31;

%z
coord3(3,1:round(94.8942/T)) = 1.03;
coord3(3,round(94.8942/T)+1:length(AA3) +round(94.8942/T)) = PP3 + 1.03;
coord3(3,length(AA3) +round(94.8942/T)+1:end) = PP3(end) + 1.03;

% position of readers
%     x1 = [-0.05, 1.5];x2 = [2, 3];x3 = [2.7, 0.05];
x1 = [0,    0,  0.756];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.07];%[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52,  0.756];%[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.07];
% x4 = [-0.3, 0.9, 0.8];

% figure
% plot3(coord3(1,:), coord3(2,:), coord3(3,:), 'LineWidth', 3);xlim([0 2.5]);ylim([0 4]);zlim([0.7 1.2])
% hold on;
% scatter3(x1(1), x1(2), x1(3), '*'); text(x1(1), x1(2) + 0.1, x1(3),"Reader 1 (0,   0,    0.705)",'fontsize',16);
% hold on;
% scatter3(x2(1), x2(2), x2(3), '*'); text(x2(1), x2(2) + 0.1, x2(3),"Reader 2 (3.2, 0,    1.105)",'fontsize',16);
% hold on;
% scatter3(x3(1), x3(2), x3(3), '*'); text(x3(1), x3(2) + 0.1, x3(3),"Reader 3 (3.2, 3.4,  0.705)",'fontsize',16);
% hold on;
% scatter3(x4(1), x4(2), x4(3), '*'); text(x4(1), x4(2) + 0.1, x4(3),"Reader 4 (0,   3.4,  1.105)",'fontsize',16);
% grid on;
% xlabel('x');ylabel('y');zlabel('z');

%save('coord3.mat','coord3')

%%
% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = 4;         % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.00001; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

H1 = NaN(1,length(coord3));     phi1 = NaN(1,length(coord3));       phi_mu1 = NaN(1,length(coord3));
H2 = NaN(1,length(coord3));     phi2 = NaN(1,length(coord3));       phi_mu2 = NaN(1,length(coord3));
H3 = NaN(1,length(coord3));     phi3 = NaN(1,length(coord3));       phi_mu3 = NaN(1,length(coord3));
H4 = NaN(1,length(coord3));     phi4 = NaN(1,length(coord3));       phi_mu4 = NaN(1,length(coord3));

r_sim1 = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3));   diff1 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3));   diff2 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3));   diff3 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3));   diff4 = NaN(1,length(coord3));

for k = 1:1:length(coord3)-1  
[H1(k+1),phi_mu1(k+1),r_sim1(k+1),rdot_sim1(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T);
[H2(k+1),phi_mu2(k+1),r_sim2(k+1),rdot_sim2(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T);
[H3(k+1),phi_mu3(k+1),r_sim3(k+1),rdot_sim3(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T);   
[H4(k+1),phi_mu4(k+1),r_sim4(k+1),rdot_sim4(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,T);  
end

