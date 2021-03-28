function J = getJacoN(y,x)

dx = 1e-12;
I = eye(length(x));
J = [ ];


for j = 1 : length(x)
    e = I(:,j);
    E = getEPA(y,x-dx*e);
    Edx = getEPA(y,x+dx*e);
    J(:,j) = (E-Edx)/(2*dx);
end
 
end