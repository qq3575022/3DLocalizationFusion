function H = getHxkPVA(x)

H = NaN(9,1);

H(1,1) = x(1)*cos(x(8))*cos(x(6))+x(2)*(cos(x(8))*sin(x(6))*sin(x(4)) - sin(x(8))*cos(x(4))) + x(3)*(cos(x(8))*sin(x(6))*cos(x(4))+sin(x(8))*sin(x(4)));
H(2,1) = x(1)*sin(x(8))*cos(x(6))+x(2)*(sin(x(8))*sin(x(6))*sin(x(4)) + cos(x(8))*cos(x(4))) + x(3)*(sin(x(8))*sin(x(6))*cos(x(4))-cos(x(8))*sin(x(4)));
H(3,1) = -x(1)*sin(x(6)) + x(2)*cos(x(6))*sin(x(4)) + x(3)*cos(x(6))*cos(x(4));

H(4,1) = x(4);
H(5,1) = x(5);
H(6,1) = x(6);
H(7,1) = x(7);
H(8,1) = x(8);
H(9,1) = x(9);
    
end