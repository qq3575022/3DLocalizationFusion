function F = getqF(gyro, Tgyro, m)
gyrox = gyro(1,m) + 0.11957;
gyroy = gyro(2,m) + 0.129038;
gyroz = gyro(3,m) - 0.114598;

T = Tgyro(m);

I   =  [1, -0.5*gyrox*T, -0.5*gyroy*T, -0.5*gyroz*T;
        0.5*gyrox*T,  1,  0.5*gyroz*T, -0.5*gyroy*T;
        0.5*gyroy*T, -0.5*gyroz*T,  1,  0.5*gyrox*T;
        0.5*gyroz*T,  0.5*gyroy*T, -0.5*gyrox*T,  1];
    
F = I;

end