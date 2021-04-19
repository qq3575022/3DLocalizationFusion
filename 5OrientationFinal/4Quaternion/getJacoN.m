function J = getJacoN(y,x,index)

dx = 1e-12;
I = eye(length(x));
J = [ ];

if index == 1

    for j = 1 : length(x)
        e = I(:,j);
        E = getEPA(y,x-dx*e);
        Edx = getEPA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
% GYRO ACC    
elseif index == 2
    
    for j = 1 : length(x)
        e = I(:,j);
        E = getEA(y,x-dx*e);
        Edx = getEA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
% MAG ACC     
elseif index == 3
    
    for j = 1 : length(x)
        e = I(:,j);
        E = getEP(y,x-dx*e);
        Edx = getEP(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
end