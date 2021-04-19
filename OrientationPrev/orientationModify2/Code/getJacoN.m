function J = getJacoN(y,x,a)

dx = 1e-12;
I = eye(length(x));
J = [ ];


for j = 1 : length(x)
    e = I(:,j);
    E = getEA(y,x-dx*e,a);
    Edx = getEA(y,x+dx*e,a);
    J(:,j) = (E-Edx)/(2*dx);
end
    


end