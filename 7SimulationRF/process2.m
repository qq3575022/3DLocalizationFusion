clc, clear, close all
% Reader 1
load('r_sim1.mat'),load('r_sim2.mat'),load('r_sim3.mat'),load('r_sim4.mat')
load('td3.mat')
figure
subplot(311), plot(r_sim1)
subplot(312), plot(r_sim2)
subplot(313), plot(r_sim3)
%%
%D1 = read_complex_binary ('/Users/Joanna/Documents/MATLAB/3D/MeasurementOct/1005/1005data/1005_reader1_xyz_4.bin', 1000000000, 1);

% Reader 1
D0_1 = read_complex_binary ('../Data/0304Reader/0304_reader1_3.bin', 1000000000, 1);
D1 = D0_1(55000:475000,1);

magD1 = abs(D1).^(-0.5);
phaseD1 = angle(D1);

time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);

%[time2_1, magD3_1, magD4_1, magD5_1, magD6_1, magD7_1, phaseD2_1] = filterAvg(magD1, phaseD1);
magD1f = medfilt1(magD1,300);
phaseD1f = medfilt1(phaseD1,300);


% Reader 2
D0_2 = read_complex_binary ('../Data/0304Reader/0304_reader2_3.bin', 1000000000, 1);
D2 = D0_2(92000:512000,1);

magD2 = abs(D2).^(-0.5);
phaseD2 = angle(D2);

time2 = 0:length(magD2)/(length(magD2)*6000):length(magD2)/6000 - length(magD2)/(length(magD2)*6000);

%[time2_2, magD3_2, magD4_2, magD5_2, magD6_2, magD7_2, phaseD2_2] = filterAvg(magD2, phaseD2);
magD2f = medfilt1(magD2,300);
phaseD2f = medfilt1(phaseD2,300);

% Reader 3
D0_3 = read_complex_binary ('../Data/0304Reader/0304_reader3_3.bin', 1000000000, 1);
D3 = D0_3(414290:885336,1);

magD3 = abs(D3).^(-0.5);
phaseD3 = angle(D3);

time3 = 0:length(magD3)/(length(magD3)*6729):length(magD3)/6729 - length(magD3)/(length(magD3)*6729);

%[time2_3, magD3_3, magD4_3, magD5_3, magD6_3, magD7_3, phaseD2_3] = filterAvg2(magD3, phaseD3);
magD3f = medfilt1(magD3,300);
phaseD3f = medfilt1(phaseD3,300);


% Reader 4
D0_4 = read_complex_binary ('../Data/0304Reader/0304_reader4_3.bin', 1000000000, 1);
D4 = D0_4(90000:510000,1);

magD4 = abs(D4).^(-0.5);
phaseD4 = angle(D4);

time4 = 0:length(magD4)/(length(magD4)*6000):length(magD4)/6000 - length(magD4)/(length(magD4)*6000);


%[time2_4, magD3_4, magD4_4, magD5_4, magD6_4, magD7_4, phaseD2_4] = filterAvg(magD4, phaseD4);

magD4f = medfilt1(magD4,300);
phaseD4f = medfilt1(phaseD4,300);



figure 
subplot(411), 
plot(time1, 0.0045*magD1+1.3), title('Radial Distantce r1')
hold on
plot(td3, r_sim1, 'LineWidth', 2), 
hold on
plot(time1, 0.0045*magD1f+1.3, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
ylim([1.5, 3]);ylabel('Distance [m]')
grid on 
grid minor


subplot(412), 
plot(time2, 0.0025*magD2+1.5), title('Radial Distantce r2')
hold on
plot(td3, r_sim2, 'LineWidth', 2), 
hold on
plot(time2, 0.0025*magD2f+1.5, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
ylim([1.5, 3.5]);ylabel('Distance [m]')
grid on 
grid minor

subplot(413), 
plot(time3, 0.0055*magD3+0.95), title('Radial Distantce r3')
hold on
plot(td3, r_sim3, 'LineWidth', 2), 
hold on
plot(time3, 0.0055*magD3f+0.95, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
ylim([1.1, 2]);ylabel('Distance [m]')
xlim([0, 70]);
grid on 
grid minor

subplot(414), 
plot(time4, 0.0015*magD4+1.2), title('Radial Distantce r4')
hold on
plot(td3, r_sim4, 'LineWidth', 2), 
hold on
plot(time4, 0.0015*magD4f+1.2, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthWest')
ylim([1, 2.5]);ylabel('Distance [m]')
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