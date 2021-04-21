clc, clear, close all

data=readtable('3.csv','Delimiter', ',');  g=9.7953;
% 1 time,   2 acc/gyro/mag,   3 x,   4 y,   5 z
% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));   data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));


acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;     acc_time = round(acc_time, 4);           
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;   gyro_time = round(gyro_time, 4);
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;     mag_time = round(mag_time, 4);

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x =  data_x(find_acc_data);   acc_data_y = data_y(find_acc_data);   acc_data_z = data_z(find_acc_data);
gyro_data_x = data_x(find_gyro_data); gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x =  data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

time = unique(sort([acc_time; gyro_time; mag_time]),'rows');

% 3.5  7.5   46  48.5   95.5  98
% x
xS = find(abs(acc_time-3.6)<0.002); xS = xS(1);
xE = find(abs(acc_time-7.5)<0.001); xE = xE(1);

xxS = find(abs(gyro_time-3.6)<0.002); xxS = xxS(1);
xxE = find(abs(gyro_time-7.5)<0.001); xxE = xxE(1);

xSS = find(abs(mag_time-3.6)<0.002); xSS = xSS(1);
xEE = find(abs(mag_time-7.5)<0.004); xEE = xEE(1);

xSSS = find(abs(time-3.6)<0.002); xSSS = xSSS(1);
xEEE = find(abs(time-7.5)<0.004); xEEE = xEEE(1);

% y
yS = find(abs(acc_time-46.1784)<0.002); yS = yS(1);
yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);

yyS = find(abs(gyro_time-46.1784)<0.002);  yyS = yyS(1);
yyE = find(abs(gyro_time-48.5)<0.002); yyE = yyE(1);

ySS = find(abs(mag_time-46.1784)<0.009);   ySS = ySS(1);
yEE = find(abs(mag_time-48.5)<0.009); yEE = yEE(1);

ySSS = find(abs(time-46.1784)<0.009);   ySSS = ySSS(1);
yEEE = find(abs(time-48.5)<0.009); yEEE = yEEE(1);

% z
zS = find(abs(acc_time-95.5027)<0.002); zS = zS(1);
zE = find(abs(acc_time-98)<0.002); zE = zE(1);

zzS = find(abs(gyro_time-95.5027)<0.002); zzS = zzS(1);
zzE = find(abs(gyro_time-98)<0.002); zzE = zzE(1);

zSS = find(abs(mag_time-95.5027)<0.009); zSS = zSS(1);
zEE = find(abs(mag_time-98)<0.009); zEE = zEE(1);

zSSS = find(abs(time-95.5027)<0.009); zSSS = zSSS(1);
zEEE = find(abs(time-98)<0.009); zEEE = zEEE(1);

accTx = acc_time(xS:xE); gyroTx = gyro_time(xxS:xxE);  magTx = mag_time(xSS:xEE);  timex = unique(sort([accTx; gyroTx; magTx]),'rows');
accTy = acc_time(yS:yE); gyroTy = gyro_time(yyS:yyE);  magTy = mag_time(ySS:yEE);  timey = unique(sort([accTy; gyroTy; magTy]),'rows');
accTz = acc_time(zS:zE); gyroTz = gyro_time(zzS:zzE);  magTz = mag_time(zSS:zEE);  timez = unique(sort([accTz; gyroTz; magTz]),'rows');



%
tdx =  timex - timex(1); %acc_time(xS:xE) - acc_time(xS);
tdy =  timey - timey(1);%acc_time(yS:yE) - acc_time(yS);
tdz =  timez - timez(1);%acc_time(zS:zE) - acc_time(zS);
%
[PP1, VV1, AA1] = groundtruth1Dx(accTx-accTx(1));
[PP2, VV2, AA2] = groundtruth1Dy(accTy-accTy(1));
[PP3, VV3, AA3] = groundtruth1Dz(accTz-accTz(1));

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

%
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

accx = medfilt1(acc_data_x,50);
accy = medfilt1(acc_data_y,50) - g;
accz = medfilt1(acc_data_z,50);

% xS = find(abs(acc_time-5)<0.002); xS = xS(1);
% BEY =find(abs(acc_time-90.5178)<0.002);BEY = BEY(1);

A3 = A3 + [0.14273;0.650848;0]*ones(1, length(acc_time));
A3_sim = A3 + [0;-0.0005;0]*acc_time';

A3_sim = awgn(A3_sim,35,'measured');

% [mean(A3_sim(1,:) - A3(1,:)), rms(A3_sim(1,:) - A3(1,:)), sqrt(rms(A3_sim(1,:) - A3(1,:)))]
% [mean(A3_sim(2,:) - A3(2,:)), rms(A3_sim(2,:) - A3(2,:)), sqrt(rms(A3_sim(2,:) - A3(2,:)))]
% [mean(A3_sim(3,:) - A3(3,:)), rms(A3_sim(3,:) - A3(3,:)), sqrt(rms(A3_sim(3,:) - A3(3,:)))]
% 
% [mean(accx' - A3(1,:)), rms(accx' - A3(1,:)), sqrt(rms(accx' - A3(1,:)))]
% [mean(accz' + A3(2,:)), rms(accz' + A3(2,:)), sqrt(rms(accz' + A3(2,:)))]
% [mean(accy' - A3(3,:)), rms(accy' - A3(3,:)), sqrt(rms(accy' - A3(3,:)))]

% figure
% subplot(611),plot(acc_time, A3_sim(1,:),'LineWidth',2), xlim([0,100]);%hold on, plot(acc_time, accx,'LineWidth',2), xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along x Axe Based on Error Model')
% subplot(612),plot(acc_time, accx,'LineWidth',2),xlim([0,100]);%hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
% subplot(613),plot(acc_time, A3_sim(2,:),'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
% subplot(614),plot(acc_time, accz,'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
% subplot(615),plot(acc_time, A3_sim(3,:),'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')
% subplot(616),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',2),xlim([0,100]); %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')

figure
subplot(311),plot(acc_time, accx,'LineWidth',3),hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(312),plot(acc_time, accz,'LineWidth',3),hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(313),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',3),hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')

%
[PP1, VV1, AA1] = groundtruth1DxRF(tdx);
[PP2, VV2, AA2] = groundtruth1DyRF(tdy);
[PP3, VV3, AA3] = groundtruth1DzRF(tdz);

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

%
% 1D - 3D coordinates
coord3 = zeros(3, length(time));
%%
% x
coord3(1,1:xSSS) = 1.03;
coord3(1,xSSS+1:length(AA1)+xSSS) = PP1 + 1.03;
coord3(1,length(AA1)+xSSS+1:end) = PP1(end) + 1.03;

%y
coord3(2,1:ySSS) = 1.31;
coord3(2,ySSS+1 : length(AA2)+ySSS) = PP2 + 1.31;
coord3(2,length(AA2)+ySSS+1:end)  = PP2(end) + 1.31;
%%
%z
coord3(3,1:zSSS) = 1.03;
coord3(3,zSSS+1 : length(AA3)+zSSS) = PP3 + 1.03;
%%
coord3(3,length(AA3)+zSSS+1:end)  = PP3(end) + 1.03;

figure
subplot(311), plot(time, coord3(1,:))
subplot(312), plot(time, coord3(2,:))
subplot(313), plot(time, coord3(3,:))
%%
% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = 4;         % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

x1 = [0,    0,  0.756];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.07];%[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52,  0.756];%[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.07];

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.00012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

zdiff1 = z(1,:)-z_prev(1,:);
zdiff2 = z(2,:)-z_prev(2,:);
zdiff3 = z(3,:)-z_prev(3,:);

%%
H1 = NaN(1,length(coord3));     phi1 = NaN(1,length(coord3));       phi_mu1 = NaN(1,length(coord3));
H2 = NaN(1,length(coord3));     phi2 = NaN(1,length(coord3));       phi_mu2 = NaN(1,length(coord3));
H3 = NaN(1,length(coord3));     phi3 = NaN(1,length(coord3));       phi_mu3 = NaN(1,length(coord3));
H4 = NaN(1,length(coord3));     phi4 = NaN(1,length(coord3));       phi_mu4 = NaN(1,length(coord3));

r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  diff1 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  diff2 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  diff3 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3)); diff4 = NaN(1,length(coord3));

for k = 1:1:length(time)-1  
[H1(k+1),phi_mu1(k+1),r_sim1(k+1),r_sim1_(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k));
[H2(k+1),phi_mu2(k+1),r_sim2(k+1),r_sim2_(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k));
[H3(k+1),phi_mu3(k+1),r_sim3(k+1),r_sim3_(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k));   
[H4(k+1),phi_mu4(k+1),r_sim4(k+1),r_sim4_(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k));  
end
%[v,         phi_mod, r,          r2,            rdot,         rdot2,           diff]

%%

%%
rdot_sim1_ = rdot_sim1 + 0.005*rand(1, length(rdot_sim1));% - 0.0025;
rdot_sim2_ = rdot_sim2 + 0.005*rand(1, length(rdot_sim2));% - 0.0025;
rdot_sim3_ = rdot_sim3 + 0.005*rand(1, length(rdot_sim3));% - 0.0025;
rdot_sim4_ = rdot_sim4 + 0.005*rand(1, length(rdot_sim4));% - 0.0025;

figure
subplot(4,1,1),plot(time, H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,2),plot(time, H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,3),plot(time, H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,4),plot(time, H4,'LineWidth',2);title('Simulated H from Reader $\#4$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')

%r_sim1_ = r_sim1+0.1*rand(1, length(r_sim1))-0.5*0.1; r_sim2_ = r_sim2+rand(1, length(r_sim2)); r_sim3_ = r_sim3+rand(1, length(r_sim3)); r_sim4_ = r_sim4+rand(1, length(r_sim4));

figure
subplot(8,1,1),plot(time, r_sim1_,'LineWidth',2);hold on;plot(time, r_sim1,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_1$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,2),plot(time, rdot_sim1_,'-','LineWidth',2);hold on;plot(time, rdot_sim1,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_1$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,3),plot(time, r_sim2_,'LineWidth',2);hold on;plot(time, r_sim2,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_2$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,4),plot(time, rdot_sim2_,'-','LineWidth',2);hold on;plot(time, rdot_sim2,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_2$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,5),plot(time, r_sim3_,'LineWidth',2);hold on;plot(time, r_sim3,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_3$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,6),plot(time, rdot_sim3_,'-','LineWidth',2);hold on;plot(time, rdot_sim3,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_3$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,7),plot(time, r_sim4_,'LineWidth',2);hold on;plot(time, r_sim4,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_4$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,8),plot(time, rdot_sim4_,'-','LineWidth',2);hold on;plot(time, rdot_sim4,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_4$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')