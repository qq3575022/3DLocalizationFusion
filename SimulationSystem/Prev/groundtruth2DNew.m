function [AXY, XX, r1, r2, r3] = groundtruth2DNew(tf, w, wf)
%The x and y coordinates of tag
%  x1 = [-0.05, 1.5];
%  x2 = [2, 3.0];
%  x3 = [2.7, 0.05];
x1 = [-1.9966,-0.0272];
x2 = [-0.0151, 1.5009];
x3 = [ 0.7439,-1.4804];
 
r = 0.3045; % radius
T = mean(diff(tf)); % sampling period

% get position x,y
XX(1,:) =  r*cos(w);
XX(2,:) = -r*sin(w).*wf;
XX(3,:) = -r*sin(w).*wf.^2 - r*cos(w).*[diff(wf)/T;NaN];

XX(4,:) =  r*sin(w);
XX(5,:) =  r*cos(w).*wf;
XX(6,:) = -r*cos(w).*wf.^2 + r*sin(w).*[diff(wf)/T;NaN];

AXY(1,:) = -r*wf.^2;
AXY(2,:) = r*[diff(wf)/T;NaN];


% get r1 r2 r3
r1 = NaN(2,length(w));  r2 = NaN(2,length(w)); r3 = NaN(2,length(w));
axx = NaN(1,length(w));ayy = NaN(1,length(w));

r1(1,1) = sqrt((XX(1,1)-x1(1))^2 + (XX(4,1)-x1(2))^2); r1(2,1) = 0;
r2(1,1) = sqrt((XX(1,1)-x2(1))^2 + (XX(4,1)-x2(2))^2); r2(2,1) = 0;
r3(1,1) = sqrt((XX(1,1)-x3(1))^2 + (XX(4,1)-x3(2))^2); r3(2,1) = 0;

for n = 2:length(tf)
    
    H = getHxk([XX(1,n),XX(2,n),XX(3,n),XX(4,n),XX(5,n),XX(6,n)]);
    r1(1,n) = H(1);
    r1(2,n) = H(2);
    r2(1,n) = H(3);
    r2(2,n) = H(4);
    r3(1,n) = H(5);
    r3(2,n) = H(6);
    axx(n)  = H(7);
    ayy(n)  = H(8);
    
end

end