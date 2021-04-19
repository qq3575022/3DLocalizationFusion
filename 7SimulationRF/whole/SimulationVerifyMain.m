clc, clear, close all

load('PropGp1DVelocityAccelerationData.mat')
ta = 0.4657;
tb = 1.1657;

td = Time(Time >= ta & Time <= tb)'-ta;
T = mean(diff(td));

td  = 0:T:13;
td1 = 0:T:23;
td3 = 0:T:13*2+23*2;

[PP1, VV1, AA1] = groundtruth1Dx(td);
[PP2, VV2, AA2] = groundtruth1Dy(td1);
[PP3, VV3, AA3] = groundtruth1Dz(td);

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

figure
subplot(311), plot(td, PP1)
subplot(312), plot(td, VV1)
subplot(313), plot(td, AA1)

figure
subplot(311), plot(td1, PP2)
subplot(312), plot(td1, VV2)
subplot(313), plot(td1, AA2)

figure
subplot(311), plot(td, PP3)
subplot(312), plot(td, VV3)
subplot(313), plot(td, AA3)
%%
% 1D - 3D coordinates
% td3
A3 =  zeros(3, 2*length(td)+2*length(td1));

% x
A3(1,1:length(AA1)) = AA1;           %      A3(1,length(AA1)*2+length(AA2)+1:length(AA1)*3+length(AA2)) = flip(AA1); 
%y
A3(2,length(AA1)+1:length(AA1) + length(AA2)) = AA2;  % A3(2,length(AA1)*2 + length(AA2)+1:end) = flip(AA2);
%z
A3(3,length(AA1) + length(AA2)+1:length(AA1)*2 + length(AA2)) = AA3;% A3(2,length(AA1)*2+length(AA2)+1:length(AA1)*3+length(AA2)) = flip(AA3);


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


% compare
data=readtable('3.csv','Delimiter', ',');  g=9.7953;

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
accx = medfilt1(acc_data_x,50);
accy = medfilt1(acc_data_y,50) - g;
accz = medfilt1(acc_data_z,50);

% xS = find(abs(acc_time-5)<0.002); xS = xS(1);
% BEY =find(abs(acc_time-90.5178)<0.002);BEY = BEY(1);
% 
% A3 = A3 + [0.14273;0.650848;0]*ones(1, 4*length(td));
% A3_sim = A3 + [0;-0.0005;0]*td3;
% 
% A3_sim = awgn(A3_sim,35,'measured');

%%

% 1D - 3D coordinates
coord3 = zeros(3, length(td3));
%%
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

coord3_ = zeros(3, length(AA2));
coord3_(1,1:length(AA1)) = flip(PP1) + 1.03;
coord3_(1,length(AA1)+1:end) = 1.03;
%y
coord3_(2,1:length(AA2)) = flip(PP2) + 1.31;

%z
coord3_(3,1:length(AA1)) = flip(PP3) + 1.03;
coord3_(3,length(AA1)+1:end) = 1.03;
%%
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
sigma = 0.00012; 

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

r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  diff1 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  diff2 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  diff3 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3)); diff4 = NaN(1,length(coord3));

for k = 1:1:length(coord3)-1  
[H1(k+1),phi_mu1(k+1),r_sim1(k+1),r_sim1_(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T);
[H2(k+1),phi_mu2(k+1),r_sim2(k+1),r_sim2_(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T);
[H3(k+1),phi_mu3(k+1),r_sim3(k+1),r_sim3_(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T);   
[H4(k+1),phi_mu4(k+1),r_sim4(k+1),r_sim4_(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,T);  
end

rdot_sim1_ = rdot_sim1 + 0.005*rand(1, length(rdot_sim1)) - 0.0025;
rdot_sim2_ = rdot_sim2 + 0.005*rand(1, length(rdot_sim2)) - 0.0025;
rdot_sim3_ = rdot_sim3 + 0.005*rand(1, length(rdot_sim3)) - 0.0025;
rdot_sim4_ = rdot_sim4 + 0.005*rand(1, length(rdot_sim4)) - 0.0025;

figure
subplot(4,1,1),plot(td3, H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,2),plot(td3, H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,3),plot(td3, H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,4),plot(td3, H4,'LineWidth',2);title('Simulated H from Reader $\#4$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')

% save('r_sim1.mat','r_sim1')
% save('r_sim2.mat','r_sim2')
% save('r_sim3.mat','r_sim3')
% save('r_sim4.mat','r_sim4')
% 
% save('td3.mat','td3')


%     figure
%     subplot(3,1,1),plot(phi_mu1,'LineWidth',2); title('Simulated Phase from Reader $\#1$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')
%     subplot(3,1,2),plot(phi_mu2,'LineWidth',2); title('Simulated Phase from Reader $\#2$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')
%     subplot(3,1,3),plot(phi_mu3,'LineWidth',2); title('Simulated Phase from Reader $\#3$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')

%r_sim1_ = r_sim1+0.1*rand(1, length(r_sim1))-0.5*0.1; r_sim2_ = r_sim2+rand(1, length(r_sim2)); r_sim3_ = r_sim3+rand(1, length(r_sim3)); r_sim4_ = r_sim4+rand(1, length(r_sim4));

figure
subplot(8,1,1),plot(td3, r_sim1_,'LineWidth',2);hold on;plot(td3, r_sim1,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_1$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,2),plot(td3, rdot_sim1_,'-','LineWidth',2);hold on;plot(td3, rdot_sim1,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_1$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,3),plot(td3, r_sim2_,'LineWidth',2);hold on;plot(td3, r_sim2,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_2$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,4),plot(td3, rdot_sim2_,'-','LineWidth',2);hold on;plot(td3, rdot_sim2,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_2$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,5),plot(td3, r_sim3_,'LineWidth',2);hold on;plot(td3, r_sim3,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_3$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,6),plot(td3, rdot_sim3_,'-','LineWidth',2);hold on;plot(td3, rdot_sim3,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_3$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,7),plot(td3, r_sim4_,'LineWidth',2);hold on;plot(td3, r_sim4,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_4$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,8),plot(td3, rdot_sim4_,'-','LineWidth',2);hold on;plot(td3, rdot_sim4,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_4$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')

%%

%     h = axes('Parent', figure());
%     hold(h, 'off');
% 
% 
%     for i = 1:1:length(PP)*4
%         scatter3(h, coord3(1,i), coord3(2,i), coord3(3,i))
%         hold on
%         plot3(coord3(1,:), coord3(2,:), coord3(3,:))
%         xlabel('x(t)')
%         ylabel('y(t)')
%         zlabel('z(t)')
%         grid on
%         pause(T/5)
%         hold on;
%         scatter3(x1(1), x1(2), x1(3), '*'); text(x1(1), x1(2) + 0.1, x1(3),"x1");
%         hold on;
%         scatter3(x2(1), x2(2), x2(3), '*'); text(x2(1), x2(2) + 0.1, x2(3),"x2");
%         hold on;
%         scatter3(x3(1), x3(2), x3(3), '*'); text(x3(1), x3(2) + 0.1, x3(3),"x3");
%         hold on;
%         scatter3(x4(1), x4(2), x4(3), '*'); text(x4(1), x4(2) + 0.1, x4(3),"x4");
%     end

