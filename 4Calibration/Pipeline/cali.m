clc, clear, close all
load('24Stpos8.mat')
load('g.mat')

% [1 x(1) x(2)][x(4) 0     0]      [x(7)]
% [0   1  x(3)][0   x(5)   0] (y + [x(8)])
% [0   0    1 ][0    0  x(6)]      [x(9)]

T = eye(3);
K = eye(3);
b = [1; 1; 1];

y = res;

x=lsqnonlin(@(xx)getEP(y, xx, g),[0.1;0.1;0.1;1;1;1;0.1;0.1;0.1]);

%save('x.mat','x')

y = [res ones(length(res),1)];

A = eye(3,4);
A(1,1) = x(4); A(1,2) = x(1)*x(5);  A(1,3) = x(2)*x(6); A(1,4) = x(7);
               A(2,2) = x(5);       A(2,3) = x(3)*x(6); A(2,4) = x(8);
                                    A(3,3) = x(6);      A(3,4) = x(9);
                                    
yy = A * y';
yR = zeros(length(y),1);

for i = 1 : 1 : length(yR)
    yR(i) = sqrt(yy(1,i)^2 + yy(2,i)^2 + yy(3,i)^2) - g;
end

% error mean
mean(yR)
% error standard deviation
sqrt(var(yR))