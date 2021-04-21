function H = getHxkPA(x)
q0 = x(1); q1 = x(2); q2 = x(3); q3 = x(4); 
B = 100; delta = -5.33/pi;
g = 9.7932;

H = NaN(6,1);

H(1,1) = (2*q1*q3-2*q0*q2)*g;
H(2,1) = (2*q0*q1+2*q2*q3)*g;
H(3,1) = (2*q0*q0+2*q3*q3-1)*g;

H(4,1) = (q0*q0+q1*q1-0.5)*B*cos(delta)+(q1*q3-q0*q2)*B*sin(delta);
H(5,1) = (q1*q2-q0*q3)*B*cos(delta)+(q0*q1+q2*q3)*B*sin(delta);
H(6,1) = (q0*q2+q1*q3)*B*cos(delta)+(q0*q0+q3*q3-0.5)*B*sin(delta);

    
end