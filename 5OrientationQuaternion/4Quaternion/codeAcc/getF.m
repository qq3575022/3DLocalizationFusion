function F = getF(gyro, Tgyro, m)
gyrox = gyro(1,m);
gyroy = gyro(2,m);
gyroz = gyro(3,m);

T = Tgyro(m);

I   =  [1, -0.5*gyrox*T, -0.5*gyroy*T, -0.5*gyroz*T;
        0.5*gyrox*T,  1,  0.5*gyroz*T, -0.5*gyroy*T;
        0.5*gyroy*T, -0.5*gyroz*T,  1,  0.5*gyrox*T;
        0.5*gyroz*T,  0.5*gyroy*T, -0.5*gyrox*T,  1];
    
F = I;

end