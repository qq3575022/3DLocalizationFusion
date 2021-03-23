clc, clear, close all

load('magT.mat');
tdmag  = round(magT,  5);

load('magx.mat'); load('magy.mat'); load('magz.mat');
mag = NaN(3,length(tdmag)); 
mag(1,:) = magx;   mag(2,:) = magy;   mag(3,:) = magz; 

magS = zeros(1,length(tdmag));
for i = 1:1:length(tdmag)
    magS(i) = mag(1,i)^2+mag(2,i)^2+mag(3,i)^2;
end
figure
plot(tdmag, magS)


%
%(mag(1) - x(1))^2 + (mag(2) - x(2))^2 + (mag(3) - x(3))^2 = B^2 = 10^4;

[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=lsqnonlin(@(xx)getEP(mag, xx),[0.1;0.1;0.1]);

mag_bias = [x(1), x(2), x(3)]