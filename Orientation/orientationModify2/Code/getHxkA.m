function H = getHxkA(x,a)
%x = [1; 0.00001; 0.00001; 0; 1; 0.00001; 0; 0; 1; 0.01; -0.01; 0.01]'
H = NaN(3,1);

H(1,1) = x(1)*(a(1) + x(7)) + x(2)*(a(2) + x(8)) + x(3)*(a(3) + x(9));
H(2,1) =                      x(4)*(a(2) + x(8)) + x(5)*(a(3) + x(9));
H(3,1) =                                           x(6)*(a(3) + x(9));

% H(1,1) = x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(14))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
% H(2,1) = x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(14))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
% H(3,1) = -x(3)*sin(x(12)) + x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));
%      
end