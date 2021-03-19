clc, clear, close all
% Reader 1
load('H1.mat'),load('H2.mat'),load('H3.mat'),load('H4.mat')
load('td3.mat')
%D1 = read_complex_binary ('/Users/Joanna/Documents/MATLAB/3D/MeasurementOct/1005/1005data/1005_reader1_xyz_4.bin', 1000000000, 1);

% Reader 1
D0_1 = read_complex_binary ('../Data/0304Reader/0304_reader1_3.bin', 1000000000, 1);
D1 = D0_1(85000:445000,1);

magD1 = abs(D1);
phaseD1 = angle(D1);

time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);

[time2_1, magD3_1, magD4_1, magD5_1, magD6_1, magD7_1, phaseD2_1] = filterAvg(magD1, phaseD1);

% figure
% subplot(311), plot(time1, magD1, 'LineWidth',2), ylim([0, 2*10^(-4)]),title('Raw Measurement Data')
% hold on, plot(time2_1, magD4_1, 'LineWidth',2), ylim([0, 2*10^(-4)]),title('Filtered Data')
% hold on,  plot(td3, 0.001*8.22*H1, 'LineWidth',2), ylim([0, 2*10^(-4)]),title('Simulated Ground Truth')
% grid on 
% grid minor
% 
% subplot(312), plot(time2_1, magD4_1, 'LineWidth',2), ylim([0, 2*10^(-4)]),title('Filtered Data')
% grid on 
% grid minor
% 
% subplot(313), plot(td3, 0.001*8.22*H1, 'LineWidth',2), ylim([0, 2*10^(-4)]),title('Simulated Ground Truth')
% grid on 
% grid minor

% Reader 2
D0_2 = read_complex_binary ('../Data/0304Reader/0304_reader2_3.bin', 1000000000, 1);
D2 = D0_2(152000:512000,1);

magD2 = abs(D2);
phaseD2 = angle(D2);

time2 = 0:length(magD2)/(length(magD2)*6000):length(magD2)/6000 - length(magD2)/(length(magD2)*6000);

[time2_2, magD3_2, magD4_2, magD5_2, magD6_2, magD7_2, phaseD2_2] = filterAvg(magD2, phaseD2);

% figure
% subplot(211), plot(time2, magD2)
% subplot(212), plot(td3, 0.001*10.61344941*H2)
% grid on 
% grid minor

% Reader 3
D0_3 = read_complex_binary ('../Data/0304Reader/0304_reader3_3.bin', 1000000000, 1);
D3 = D0_3(200000:920000,1);

magD3 = abs(D3);
phaseD3 = angle(D3);

% figure
% subplot(211), plot(magD3)
% subplot(212), plot(H3)
% grid on 
% grid minor

time3 = 0:length(magD3)/(length(magD3)*12000):length(magD3)/12000 - length(magD3)/(length(magD3)*12000);


[time2_3, magD3_3, magD4_3, magD5_3, magD6_3, magD7_3, phaseD2_3] = filterAvg2(magD3, phaseD3);

% Reader 4
D0_4 = read_complex_binary ('../Data/0304Reader/0304_reader4_3.bin', 1000000000, 1);
D4 = D0_4(150000:510000,1);

magD4 = abs(D4);
phaseD4 = angle(D4);

% figure
% subplot(211), plot(magD4)
% subplot(212), plot(H4)
% grid on 
% grid minor
%%
time4 = 0:length(magD4)/(length(magD4)*6000):length(magD4)/6000 - length(magD4)/(length(magD4)*6000);


[time2_4, magD3_4, magD4_4, magD5_4, magD6_4, magD7_4, phaseD2_4] = filterAvg(magD4, phaseD4);


figure 
subplot(411), 
plot(time1, magD1), title('Reader 1')
hold on
plot(td3, 0.001*8.22*H1, 'LineWidth', 2), 
hold on
plot(time2_1, magD3_1, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
grid on 
grid minor


subplot(412), 
plot(time2, magD2), title('Reader 2')
hold on
plot(td3, 0.001*10.61344941*H2, 'LineWidth', 2), 
hold on
plot(time2_2, magD4_2, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
grid on 
grid minor

subplot(413), 
plot(time3, magD3), title('Reader 3')
hold on
plot(td3, 0.5*(0.001*42.019839*H3 - 3.70345*10^(-4)), 'LineWidth', 2), 
hold on
plot(time2_3, magD4_3, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
grid on 
grid minor

subplot(414), 
plot(time4, magD4), title('Reader 4')
hold on
plot(td3, 0.001*3.779274583*H3, 'LineWidth', 2), 
hold on
plot(time2_4, magD4_4, 'LineWidth', 2),
legend('Whole','Measurement','Simulation','Filtered','location','NorthWest')
grid on 
grid minor

% subplot(413), plot(phaseD1), title('Raw Phase Measurement')
% grid on 
% grid minor
% 
% subplot(414), plot(phaseD2), title('Averaged Phase Measurement')
% grid on 
% grid minor

%
% mean(d4(1224:1924))