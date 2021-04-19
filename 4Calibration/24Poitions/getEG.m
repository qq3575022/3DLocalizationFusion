function E = getEG(yy, xx, res, resgT, gyro_x, gyro_y, gyro_z, windows)

angles = zeros(length(windows)/2 - 1, 3);

for j = 2 : 2 : length(windows) - 2

for i = windows(j) : 1 : windows(j+1)
    [gyroX, gyroY, gyroZ] = caliGyro(gyro_x(i), gyro_y(i), gyro_z(i), xx);
    
    angles(j/2,1) = angles(j/2,1) + gyroX * (resgT(i + 1) - resgT(i));
    angles(j/2,2) = angles(j/2,2) + gyroY * (resgT(i + 1) - resgT(i));
    angles(j/2,3) = angles(j/2,3) + gyroZ * (resgT(i + 1) - resgT(i));
end

end


trans = zeros(length(angles),3);

for l = 1 : 1 : length(res) - 1
x(1) = res(l,1);    x(2) = res(l,2);    x(3) = res(l,3);
x(4) = angles(l,1); x(6) = angles(l,2); x(8) = angles(l,3);

trans(l,1) = x(1)*cos(x(8))*cos(x(6))+x(2)*(cos(x(8))*sin(x(6))*sin(x(4)) - sin(x(8))*cos(x(4))) + x(3)*(cos(x(8))*sin(x(6))*cos(x(4))+sin(x(8))*sin(x(4)));
trans(l,2) = x(1)*sin(x(8))*cos(x(6))+x(2)*(sin(x(8))*sin(x(6))*sin(x(4)) + cos(x(8))*cos(x(4))) + x(3)*(sin(x(8))*sin(x(6))*cos(x(4))-cos(x(8))*sin(x(4)));
trans(l,3) = -x(1)*sin(x(6)) + x(2)*cos(x(6))*sin(x(4)) + x(3)*cos(x(6))*cos(x(4));
end

% 
E = 0;

for j = 1 : 1 : length(trans)
    E = E + (trans(j,1) - yy(1,j+1))^2 + (trans(j,2) - yy(2,j+1))^2 + (trans(j,3) - yy(3,j+1))^2;
end

end