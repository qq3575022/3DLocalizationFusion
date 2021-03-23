function H = getHxkA(x)

H = NaN(3,1);

H(1,1) = x(1)*cos(x(8))*cos(x(6))+x(2)*(cos(x(8))*sin(x(6))*sin(x(4)) - sin(x(8))*cos(x(4))) + x(3)*(cos(x(8))*sin(x(6))*cos(x(4))+sin(x(8))*sin(x(4)));
H(2,1) = x(1)*sin(x(8))*cos(x(6))+x(2)*(sin(x(8))*sin(x(6))*sin(x(4)) + cos(x(8))*cos(x(4))) + x(3)*(sin(x(8))*sin(x(6))*cos(x(4))-cos(x(8))*sin(x(4)));
H(3,1) = -x(1)*sin(x(6)) + x(2)*cos(x(6))*sin(x(4)) + x(3)*cos(x(6))*cos(x(4));
    
end