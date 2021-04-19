function J = getJacoN2(x)

J = zeros(4);
J(1,3) = 1;
J(2,4) = 1; 
J(3,1) =  cos(x(3));
J(3,2) = -sin(x(3));
J(3,3) = -x(1)*sin(x(3)) - x(2)*cos(x(3));

J(4,1) =  sin(x(3));
J(4,2) =  cos(x(3));
J(4,3) =  x(1)*cos(x(3)) - x(2)*sin(x(3));

end