function J = getJacoN(y,x,index)

dx = 1e-12;
I = eye(length(x));
J = [ ];

if index == 1
    for j = 1 : length(x)
        e = I(:,j);
        E = getEMGA(y,x-dx*e);
        Edx = getEMGA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end

elseif index == 2
    for j = 1 : length(x)
        e = I(:,j);
        E = getEGA(y,x-dx*e);
        Edx = getEGA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 3
    for j = 1 : length(x)
        e = I(:,j);
        E = getEMG(y,x-dx*e);
        Edx = getEMG(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 4
    for j = 1 : length(x)
        e = I(:,j);
        E = getEMA(y,x-dx*e);
        Edx = getEMA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 5
    for j = 1 : length(x)
        e = I(:,j);
        E = getEA(y,x-dx*e);
        Edx = getEA(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 6
    for j = 1 : length(x)
        e = I(:,j);
        E = getEG(y,x-dx*e);
        Edx = getEG(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end

elseif index == 7
    for j = 1 : length(x)
        e = I(:,j);
        E = getEM(y,x-dx*e);
        Edx = getEM(y,x+dx*e);
        J(:,j) = (E-Edx)/(2*dx);
    end
    

end