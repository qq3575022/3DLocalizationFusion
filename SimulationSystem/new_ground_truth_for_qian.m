clc, clear, close all

%% experimental data

load('Measurement/yVectorData.mat')
%%
t = tf;
T = mean(diff(t));

r1 = r1f;
r2 = r2f;
r3 = r3f;  

r1dot = r1dot_e;
r2dot = r2dot_e;
r3dot = r3dot_e;

psif1 = unwrap(psif1)-2*pi;
theta = psif1+3*pi/2;
thetadot = wf;

ax = ax1+0.62;
ay = ay1+0.26;

figure
plot(t,wf)

%% ground truth

rho = 0.3045;
x1 = -1.9966;
y1 = -0.0272;
x2 = -0.0151;
y2 = 1.5009;
x3 = 0.7439;
y3 = -1.4804;

x = rho*cos(theta);
y = rho*sin(theta);

xdot = -y.*thetadot;
ydot = x.*thetadot;

r1_t = sqrt((x-x1).^2+(y-y1).^2);
r2_t = sqrt((x-x2).^2+(y-y2).^2);
r3_t = sqrt((x-x3).^2+(y-y3).^2);

r1dot_t = ((x-x1).*xdot+(y-y1).*ydot)./r1_t;
r2dot_t = ((x-x2).*xdot+(y-y2).*ydot)./r2_t;
r3dot_t = ((x-x3).*xdot+(y-y3).*ydot)./r3_t;

ax_t = -rho*thetadot.^2;
ay_t = rho*[diff(thetadot)/T;NaN];

 errormea1 = r1-r1_t;   errormea2 = r2-r2_t;   errormea3 = r3-r3_t;
   
error_matrics_meas1 = [mean(errormea1), var(errormea1),sqrt(var(errormea1))]
error_matrics_meas2 = [mean(errormea2), var(errormea2),sqrt(var(errormea2))]
error_matrics_meas3 = [mean(errormea3), var(errormea3),sqrt(var(errormea3))]
    
figure
subplot(421)
plot(t,r1,t,r1_t)
subplot(423)
plot(t,r2,t,r2_t)
subplot(425)
plot(t,r3,t,r3_t)
subplot(422)
plot(t,r1dot,t,r1dot_t)
subplot(424)
plot(t,r2dot,t,r2dot_t)
subplot(426)
plot(t,r3dot,t,r3dot_t)
subplot(427)
plot(t,ax,t,ax_t)
subplot(428)
plot(t,ay,t,ay_t)