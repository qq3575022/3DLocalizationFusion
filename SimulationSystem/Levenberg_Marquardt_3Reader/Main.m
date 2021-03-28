clc, clear, close all
sim = 1;

if sim == 1
load('PropGp1DPositionData.mat')
load('PropGp1DVelocityAccelerationData.mat')

ta = 0.4657;
tb = 1.1657;

% Can be used as measurement
pp = Position_1D(Time >= ta & Time <= tb)';
vv = MeasVel1D(Time >= ta & Time <= tb)';
aa = MeasAccel1D(Time >= ta & Time <= tb)';

td = Time(Time >= ta & Time <= tb)'-ta;
global T;
T = mean(diff(td));

[PP, VV, AA] = groundtruth1D(td);
% Start from 0
PP = PP - PP(1);

% 1D - 3D coordinates
coord3 = zeros(3, length(PP)*4);
% x
coord3(1,1:length(PP)) = PP;
coord3(1,length(PP):length(PP)*3) = PP(end);
coord3(1,length(PP)*3+1:length(PP)*4) = flip(PP);
%y
coord3(2,length(PP)+1:length(PP)*2) = PP;
coord3(2,length(PP)*2:length(PP)*3) = PP(end);
coord3(2,length(PP)*3+1:length(PP)*4) = flip(PP);
%z
coord3(3,length(PP)*2+1:length(PP)*3) = PP;
coord3(3,length(PP)*3+1:length(PP)*4) = flip(PP);

tf = NaN(1,length(coord3));
tf(1:length(td)) = td;
tf(length(td) + 1:length(td)*2) = td + tf(length(td));
tf(length(td)*2+1:length(td)*3) = td + tf(length(td)*2);
tf(length(td)*3+1:length(td)*4) = td + tf(length(td)*3);

% Acceleration
acc3 = zeros(3, length(AA)*4);

% x
acc3(1,1:length(AA)) = AA;
acc3(1,length(AA)*3+1:length(AA)*4) = flip(AA);
%y
acc3(2,length(AA)+1:length(AA)*2) = AA;
acc3(2,length(AA)*3+1:length(AA)*4) = flip(AA);
%z
acc3(3,length(AA)*2+1:length(AA)*3) = AA;
acc3(3,length(AA)*3+1:length(AA)*4) = flip(AA);

% GyroMag
gyro3 = zeros(3, length(AA)*4);
mag3  = zeros(3, length(AA)*4);

% position of readers
%     x1 = [-0.05, 1.5];x2 = [2, 3];x3 = [2.7, 0.05];
x1 = [-0.3,-0.3, 0];
x2 = [ 0.9, 0.9, 0];
x3 = [ 0.9,-0.3,-0.3];
%x4 = [-0.3, 0.9, 0.3];

% figure
% plot3(coord3(1,:), coord3(2,:), coord3(3,:));
% hold on;
% scatter3(x1(1), x1(2), x1(3), '*'); text(x1(1), x1(2) + 0.1, x1(3),"x1");
% hold on;
% scatter3(x2(1), x2(2), x2(3), '*'); text(x2(1), x2(2) + 0.1, x2(3),"x2");
% hold on;
% scatter3(x3(1), x3(2), x3(3), '*'); text(x3(1), x3(2) + 0.1, x3(3),"x3");
% hold on;
% scatter3(x4(1), x4(2), x4(3), '*'); text(x4(1), x4(2) + 0.1, x4(3),"x4");
% grid on;
% xlabel('x');ylabel('y');zlabel('z');

% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = 4;         % load modulation factor of the tag
f1 = 5.78*10^9;
f2 = 5.79*10^9;
f3 = 5.80*10^9;
f4 = 5.81*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0002; 

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

r1 = NaN(1,length(coord3)); rdot1 = NaN(1,length(coord3));   diff1 = NaN(1,length(coord3));
r2 = NaN(1,length(coord3)); rdot2 = NaN(1,length(coord3));   diff2 = NaN(1,length(coord3));
r3 = NaN(1,length(coord3)); rdot3 = NaN(1,length(coord3));   diff3 = NaN(1,length(coord3));
r4 = NaN(1,length(coord3)); rdot4 = NaN(1,length(coord3));   diff4 = NaN(1,length(coord3));

for k = 1:1:length(coord3)-1  
 [H1(k+1),phi_mu1(k+1),r1(k),rdot1(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T);
 [H2(k+1),phi_mu2(k+1),r2(k),rdot2(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T);
 [H3(k+1),phi_mu3(k+1),r3(k),rdot3(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T);   
 %[H4(k+1),phi_mu4(k+1),r4(k),rdot4(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,T);  
end
 
% Set timestep
time_step    = 0.1;   iternum    = 200; 
time_stepPos = 0.1;   iternumPos = 200;  

% ++++++++++++++++++++++   Input: Length of NLS estimation  ++++++++++++++++
meanrr = NaN(6,40);  stdrr = NaN(6,40);  time = linspace(1,40,40);
meanx = NaN(6,40);   stdx = NaN(6,40);
meany = NaN(6,40);   stdy = NaN(6,40);

%%
figure
subplot(3,1,1), plot(r1)
subplot(3,1,2), plot(r2)
subplot(3,1,3), plot(r3)
%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%==============================  Estimation Based on R ===============================
%for N = 1:1:40
N = 7;
% Simulated instantaneous trilateration y
y = NaN(3*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1,r2,r3,l,N);
end

% Initiate x, e and temproary variable
x = [.5;0;0; .5;0;0; 3;0;0; 0;0; 0;0; 0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[.5;0;0; .5;0;0; 3;0;0; 0;0; 0;0; 0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,1,T,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,1,T,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,1,T,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,1,T,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,T,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]^(N-1)*x;

errorestx = x(1,:)-coord3(1,N+1:length(tf)); erroresty = x(4,:)-coord3(2,N+1:length(tf)); 

meanx(1,N) = mean(errorestx);      meany(1,N) = mean(erroresty);
stdx(1,N) = sqrt(var(errorestx));  stdy(1,N) =  sqrt(var(erroresty));

meanrr(1,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(1,N)  = sqrt(var(errorestx)    + var(erroresty));

% ---------------------------------------------
% Get instantaneous trilateration error

% errormeasx = x_meas(1,N+1:length(tf))-XX(1,N+1:length(tf)); 
% errormeasy = x_meas(4,N+1:length(tf))-XX(4,N+1:length(tf)); 
% 
% fprintf('------------- measurement error ------------------')
% error_matricsmeas  = [sqrt(mean(errormeasx)^2+mean(errormeasy)^2), sqrt(var(errormeasx)^2+var(errormeasy)^2),sqrt(var(errormeasx) + var(errormeasy))]
% error_matricsmeasx = [mean(errormeasx), var(errormeasx),sqrt(var(errormeasx))]
% error_matricsmeasy = [mean(errormeasy), var(errormeasy),sqrt(var(errormeasy))]


%end
%%
figure
subplot(3,1,1),
%plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),coord3(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)), x(1,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)), 1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$x[m]$','interpreter','latex');
legend('ground truth','estimated','initialization');title('Position x')

subplot(3,1,2),
%plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),coord3(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),1*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$y[m]$','interpreter','latex');
legend('ground truth','estimated','initialization');title('Position y');xlabel('t [s]')

subplot(3,1,3),
plot(tf(N+1:length(tf)),coord3(3,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$\dot y[m]$','interpreter','latex');
legend({'ground truth','estimated','initialization'},'Location','NorthEast');title('Velocity y');xlabel('t [s]')

%%%%==============================================================================================================
figure
subplot(3,1,1),
%plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
%plot(tf(N+1:length(tf)),coord3(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)), x(2,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)), 1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$x[m]$','interpreter','latex');
legend('estimated','initialization');title('Position x')

subplot(3,1,2),
%plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
%plot(tf(N+1:length(tf)),coord3(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),1*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$y[m]$','interpreter','latex');
legend('estimated','initialization');title('Position y');xlabel('t [s]')

subplot(3,1,3),
%plot(tf(N+1:length(tf)),coord3(3,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',1)
ylabel('$\dot y[m]$','interpreter','latex');
legend({'ground truth','estimated','initialization'},'Location','NorthEast');title('Velocity y');xlabel('t [s]')

% %%
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % =============================  Estimation Based on R Rdot ==========================
% for N = 1:1:40
% % Simulated instantaneous trilateration y
% y = NaN(6*N,length(tf)-N);
% 
% % Get yN
% for l = 1:1:length(tf)-N
% y(:,l) = getyNPV(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,l,N);
% end
% 
% % Initiate x, e and temproary variable
% x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);
% 
% for m = 1:1:length(tf)-N %time 1 t
%     
% x(:,m)=lsqnonlin(@(xx)getEPV(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
% end
% 
% % Get F^(N-1)*x
% x = [1,T,0,0,0,0,0,0,0;
%      0,1,T,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,0,0;
%      0,0,0,1,T,0,0,0,0;
%      0,0,0,0,1,T,0,0,0;
%      0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,1,T,0;
%      0,0,0,0,0,0,0,1,T;
%      0,0,0,0,0,0,0,0,1]^(N-1)*x;
% 
% % Get estimation error
% errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 
% 
% % error_matricsestxRRdot = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
% % error_matricsestyRRdot = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
% % fprintf('-------------------------------')
% meanx(2,N) = mean(errorestx);     meany(2,N) = mean(erroresty);
% stdx(2,N) = sqrt(var(errorestx)); stdy(2,N) =  sqrt(var(erroresty));
% 
% meanrr(2,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
% stdrr(2,N)  = sqrt(var(errorestx)    + var(erroresty));
% end
% 
% %%
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % ==========================  Estimation R Acc Orientation ===========================
% for N = 1:1:40
% y = NaN(N*6,length(tf)-N);
% 
% % Get yN
% for l = 1:1:length(tf)-N
% y(:,l) = getyNPAO(r1f,r2f,r3f,phi(1,:),AXY(1,:),AXY(2,:),l,N);
% end
% 
% % Initiate x, e and temproary variable
% x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);
% 
% for m = 1:1:length(tf)-N %time 1 t
%     
% x(:,m)=lsqnonlin(@(xx)getEPAO(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
% end
% 
% % Get F^(N-1)*x
% x = [1,T,0,0,0,0,0,0,0;
%      0,1,T,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,0,0;
%      0,0,0,1,T,0,0,0,0;
%      0,0,0,0,1,T,0,0,0;
%      0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,1,T,0;
%      0,0,0,0,0,0,0,1,T;
%      0,0,0,0,0,0,0,0,1]^(N-1)*x;
% 
% % Get estimation error
% errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 
% 
% meanx(3,N) = mean(errorestx);     meany(3,N) = mean(erroresty);
% stdx(3,N) = sqrt(var(errorestx)); stdy(3,N) =  sqrt(var(erroresty));
% 
% meanrr(3,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
% stdrr(3,N)  = sqrt(var(errorestx)    + var(erroresty));
% end
% 
% %%
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % =============================  Estimation Based on R Acc ===========================
% for N = 1:1:40
% y =   NaN(7*N,length(tf)-N); 
% 
% % Get yN
% for l = 1:1:length(tf)-N
% y(:,l) = getyNPA(r1f,r2f,r3f,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
% end
% 
% % Initiate x, e and temproary variable
% x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);
% 
% for m = 1:1:length(tf)-N %time 1 t
% x(:,m)=lsqnonlin(@(xx)getEPA(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
% end
% 
% % Get F^(N-1)*x
% x = [1,T,0,0,0,0,0,0,0;
%      0,1,T,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,0,0;
%      0,0,0,1,T,0,0,0,0;
%      0,0,0,0,1,T,0,0,0;
%      0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,1,T,0;
%      0,0,0,0,0,0,0,1,T;
%      0,0,0,0,0,0,0,0,1]^(N-1)*x;
% 
% % Get estimation error
% errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 
% 
% meanx(4,N) = mean(errorestx);     meany(4,N) = mean(erroresty);
% stdx(4,N) = sqrt(var(errorestx)); stdy(4,N) =  sqrt(var(erroresty));
% 
% meanrr(4,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
% stdrr(4,N)  = sqrt(var(errorestx)    + var(erroresty));
% 
% end
% 
% %%
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % ==========================  Estimation R Rdot Acc Orientation ======================
% for N = 1:1:40
% y = NaN(N*9,length(tf)-N);
% 
% % Get yN
% for l = 1:1:length(tf)-N
% y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
% end
% 
% % Initiate x, e and temproary variable
% x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);
% 
% for m = 1:1:length(tf)-N %time 1 t
%     
% x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
% end
% 
% % Get F^(N-1)*x
% x = [1,T,0,0,0,0,0,0,0;
%      0,1,T,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,0,0;
%      0,0,0,1,T,0,0,0,0;
%      0,0,0,0,1,T,0,0,0;
%      0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,1,T,0;
%      0,0,0,0,0,0,0,1,T;
%      0,0,0,0,0,0,0,0,1]^(N-1)*x;
% 
% % Get estimation error
% errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 
% 
% meanx(5,N) = mean(errorestx);     meany(5,N) = mean(erroresty);
% stdx(5,N) = sqrt(var(errorestx)); stdy(5,N) =  sqrt(var(erroresty));
% 
% meanrr(5,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
% stdrr(5,N)  = sqrt(var(errorestx)    + var(erroresty));
% end
% 
% 
% %%
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % ==========================  Estimation R Rdot Acc ==================================
% for N = 1:1:40
% y = NaN(N*10,length(tf)-N);
% 
% % Get yN
% for l = 1:1:length(tf)-N
% y(:,l) = getyNPVA(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
% end
% 
% % Initiate x, e and temproary variable
% x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);
% 
% for m = 1:1:length(tf)-N %time 1 t
%     
% x(:,m)=lsqnonlin(@(xx)getEPVA(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
% end
% 
% % Get F^(N-1)*x
% x = [1,T,0,0,0,0,0,0,0;
%      0,1,T,0,0,0,0,0,0;
%      0,0,1,0,0,0,0,0,0;
%      0,0,0,1,T,0,0,0,0;
%      0,0,0,0,1,T,0,0,0;
%      0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,1,T,0;
%      0,0,0,0,0,0,0,1,T;
%      0,0,0,0,0,0,0,0,1]^(N-1)*x;
% 
% % Get estimation error
% errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 
% 
% meanx(6,N) = mean(errorestx);     meany(6,N) = mean(erroresty);
% stdx(6,N) = sqrt(var(errorestx)); stdy(6,N) =  sqrt(var(erroresty));
% 
% meanrr(6,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
% stdrr(6,N)  = sqrt(var(errorestx)    + var(erroresty));
% end
% 
% %%
% errormeasx = x_meas(1,:) - XX(1,:); errormeasy = x_meas(4,:) - XX(4,:);
% 
% time = time(1:40);
% figure
% subplot(2,1,1), plot(time,meanrr(1,1:40),time,meanrr(2,1:40),time,meanrr(3,1:40),time,meanrr(4,1:40),time,meanrr(5,1:40),time,meanrr(6,1:40),'LineWidth',2); hold on; plot(time, ones(1,40)*sqrt(mean(errormeasx)^2+mean(errormeasy)^2),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% grid on; xlim([1,40]);xlabel('Stack Length');ylim([0,0.03])
% legend('estimated r','estimated r $\dot r$','estimated r $a_x$ $a_y$ $\theta_z$','estimated r $a_x$ $a_y$ $\theta_z$ $\omega_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$ $\omega_z$','interpreter','latex','Location','NorthWest');title('Mean Error of Estimation');ylabel('Mean Error [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
% 
% subplot(2,1,2), plot(time,stdrr(1,1:40),time,stdrr(2,1:40),time,stdrr(3,1:40),time,stdrr(4,1:40),time,stdrr(5,1:40),time,stdrr(6,1:40),'LineWidth',2);hold on; plot(time, ones(1,40)*sqrt(var(errormeasx) + var(errormeasy)),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% grid on; xlim([1,40]);xlabel('Stack Length');ylim([0.015,0.075]);
% legend('estimated r','estimated r $\dot r$','estimated r $a_x$ $a_y$ $\theta_z$','estimated r $a_x$ $a_y$ $\theta_z$ $\omega_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$ $\omega_z$','instantaneous trilateration','interpreter','latex','Location','NorthWest');
% title('RMS Error of Estimation');ylabel('RMS Error [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
% 
% %%
% figure
% plot(time,stdrr(1,1:40),time,stdrr(2,1:40),time,stdrr(3,1:40),time,stdrr(4,1:40),time,stdrr(5,1:40),time,stdrr(6,1:40),'LineWidth',2);hold on; plot(time, ones(1,40)*sqrt(var(errormeasx) + var(errormeasy)),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% grid on; xlim([1,40]);xlabel('Stack Length','FontSize',22);ylim([0.015,0.075]);
% legend('estimated r','estimated r $\dot r$','estimated r $a_x$ $a_y$ $\theta_z$','estimated r $a_x$ $a_y$ $\theta_z$ $\omega_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$ $\omega_z$','instantaneous trilateration','interpreter','latex','Location','NorthWest','FontSize',22,'Position',[0.2, 0.6, .25, .25]);
% title('RMS Error of Estimation','FontSize',22);ylabel('RMS Error [m]','FontSize',22);xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
% 
% %%
% 
% save('meanrr.mat','meanrr')
% save('stdrr.mat', 'stdrr')
% 
% save('meanx_3.mat','meanx')
% save('stdx_3.mat','stdx')
% 
% save('meany_3.mat','meany')
% save('stdy_3.mat','stdy')
% 
% %%
% fprintf('----------- E ----------')
% meanE = sqrt(mean(x_meas(1,:) - XX(1,:))^2 + mean(x_meas(4,:) - XX(4,:))^2)
% varE = sqrt(var(x_meas(1,:) - XX(1,:)) + var(x_meas(4,:) - XX(4,:)))
% 
% fprintf('----------- R ----------')
% meanR = meanrr(1,10)
% minxR = min(meanrr(1,:))
% meanR = find(meanrr(1,:)==minxR)
% 
% varR  = stdrr(1,10)
% minyR = min(stdrr(1,:)) 
% varR  = find(stdrr(1,:)==minyR)
% 
% fprintf('----------- R R----------')
% meanRR = meanrr(2,10)
% minxRR = min(meanrr(2,:))
% meanRR = find(meanrr(2,:)==minxRR)
% 
% varRR  = stdrr(2,10)
% minyRR = min(stdrr(2,:)) 
% varRR  = find(stdrr(2,:)==minyRR)
% 
% fprintf('----------- R ACC O----------')
% meanRAO = meanrr(3,10)
% minxRAO = min(meanrr(3,:))
% meanRAO = find(meanrr(3,:)==minxRAO)
% 
% varRAO  = stdrr(3,10)
% minyRAO = min(stdrr(3,:)) 
% varRAO  = find(stdrr(3,:)==minyRAO)
% 
% fprintf('----------- R ACC ----------')
% meanRA = meanrr(4,10)
% minxRA = min(meanrr(4,:))
% meanRA = find(meanrr(4,:)==minxRA)
% 
% varRA  = stdrr(4,10)
% minyRA = min(stdrr(4,:)) 
% varRA  = find(stdrr(4,:)==minyRA)
% 
% fprintf('----------- R R ACC O----------')
% meanRRAO = meanrr(5,10)
% minxRRAO = min(meanrr(5,:))
% meanRRAO = find(meanrr(5,:)==minxRRAO)
% 
% varRRAO  = stdrr(5,10)
% minyRRAO = min(stdrr(5,:))  
% varRRAO  = find(stdrr(5,:)==minyRRAO)
% 
% fprintf('----------- R R ACC ----------')
% meanRRA = meanrr(6,10)
% minxRRA = min(meanrr(6,:)) 
% meanRRA = find(meanrr(6,:)==minxRRA)
% 
% varRRA  = stdrr(6,10)
% minyRRA = min(stdrr(6,:))  
% varRRA  = find(stdrr(6,:)==minyRRA)
% 

else
    
% -----------Load Data Infromation------------
% time: tf    
% radical distances: r1f, r2f, r3f
% radical velocity:  r1dot_e, r2dot_e, r3dot_e
% orientation    psif1
% angular velocity  wf
% acceleration ax1 ay1
% ---------------------------------------------
% Ground truth coordinates of the moving tag path

% ---------------------------------------------
% Get acceleration along x^B y^B

% ---------------------------------------------
% Get ground-truth

    load('yVectorData.mat')

    offset = 0;
    for n = 2:1:length(psif1)
        if psif1(n) - psif1(n-1) < -2
            offset = 2*pi;
        end
        psif1(n) = psif1(n) + offset - pi;  
    end

    psif1(1) = psif1(1) - pi;


    h = 1e-5;           tc = 0:h:2.00000;
    T = mean(diff(tf)); td = 0:T:2.00000;

    [W_A, W_V, W, XX, rr] = groundtruth2D(tc, td);

    figure
    subplot(6,1,1), plot(XX(1,:));
    subplot(6,1,2), plot(XX(2,:));
    subplot(6,1,3), plot(XX(3,:));
    subplot(6,1,4), plot(XX(4,:));
    subplot(6,1,5), plot(XX(5,:));
    subplot(6,1,6), plot(XX(6,:));
    
    
    coord3 = zeros(3, 282);
    coord3(1,1:142) = XX(1,1:142);
    coord3(1,142+1:282) = XX(1,142);
    coord3(2,:) = XX(4,1:282);
    coord3(3,142+1:282) = XX(1,142+1:282)-XX(1,142);
    
    tf = td;

    x1 = [-0.3,-0.3, 0];
    x2 = [ 0.9, 0.9, 0];
    x3 = [ 0.9,-0.3,-0.3];
    x4 = [-0.3, 0.9, 0.3];
    
    figure
    plot3(coord3(1,:), coord3(2,:), coord3(3,:));
    hold on;
    scatter3(x1(1), x1(2), x1(3), '*'); text(x1(1), x1(2) + 0.1, x1(3),"x1");
    hold on;
    scatter3(x2(1), x2(2), x2(3), '*'); text(x2(1), x2(2) + 0.1, x2(3),"x2");
    hold on;
    scatter3(x3(1), x3(2), x3(3), '*'); text(x3(1), x3(2) + 0.1, x3(3),"x3");
    hold on;
    scatter3(x4(1), x4(2), x4(3), '*'); text(x4(1), x4(2) + 0.1, x4(3),"x4");
    grid on;
    xlabel('x');ylabel('y');zlabel('z');

    %ax_noisy = AXY(1,:) + 0.2*randn(size(AXY(1,:)));ay_noisy = AXY(2,:) + 0.2*randn(size(AXY(2,:)));

    % Parameters of tag
    Gt = 14.62;    % tag's antenna gain
    X = 0.85;      % polarization mismatch
    M = 4;         % load modulation factor of the tag
    f1 = 5.78*10^9;
    f2 = 5.79*10^9;
    f3 = 5.80*10^9;
    f4 = 5.81*10^9;

    % Parameters of reader
    PT = 1;         % reader's transmitted power
    GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
    GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
    R = 15;

    % Channel noise error covariance
    sigma = 0.0002; 

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

    r1 = NaN(1,length(coord3)); rdot1 = NaN(1,length(coord3));   diff1 = NaN(1,length(coord3));
    r2 = NaN(1,length(coord3)); rdot2 = NaN(1,length(coord3));   diff2 = NaN(1,length(coord3));
    r3 = NaN(1,length(coord3)); rdot3 = NaN(1,length(coord3));   diff3 = NaN(1,length(coord3));
    r4 = NaN(1,length(coord3)); rdot4 = NaN(1,length(coord3));   diff4 = NaN(1,length(coord3));


    for k = 1:1:length(coord3)-1  
     [H1(k+1),phi_mu1(k+1),r1(k),rdot1(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T);
     [H2(k+1),phi_mu2(k+1),r2(k),rdot2(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T);
     [H3(k+1),phi_mu3(k+1),r3(k),rdot3(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T);   
     [H4(k+1),phi_mu4(k+1),r4(k),rdot4(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,T);  
    end
    
    %%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%==============================  Estimation Based on R ===============================
%for N = 1:1:40
N = 3;
% Simulated instantaneous trilateration y
y = NaN(4*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1,r2,r3,r4,l,N);
end

% Initiate x, e and temproary variable
x = [.5;0;0;0;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[0.3;0;0;0.3;0;0;0;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

errorestx = x(1,:)-coord3(1,N+1:length(tf)); erroresty = x(4,:)-coord3(2,N+1:length(tf)); 

meanx(1,N) = mean(errorestx);      meany(1,N) = mean(erroresty);
stdx(1,N) = sqrt(var(errorestx));  stdy(1,N) =  sqrt(var(erroresty));

meanrr(1,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(1,N)  = sqrt(var(errorestx)    + var(erroresty));

%end
%%
subplot(3,1,1),
%plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),coord3(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)), 1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$x[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(3,1,2),
%plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),coord3(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$y[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y');xlabel('t [s]')

subplot(3,1,3),
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),coord3(3,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$\dot y[m]$','interpreter','latex');
legend({'estimated','ground truth','initialization'},'Location','NorthEast');title('Velocity y');xlabel('t [s]')

end