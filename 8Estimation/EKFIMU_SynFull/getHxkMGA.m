function H = getHxkMGA(x)

H = NaN(9,1);


H(1,1) = x(10);
H(2,1) = x(11);
H(3,1) = x(12);
H(4,1) = x(13);
H(5,1) = x(14);
H(6,1) = x(15);

H(7,1) =  x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(14))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
H(8,1) =  x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(14))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
H(9,1) = -x(3)*sin(x(12))           +x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));

    
end