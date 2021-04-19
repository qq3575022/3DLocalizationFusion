function [gyroX, gyroY, gyroZ] = caliGyro(gyro_x, gyro_y, gyro_z, xx)
    %gyM = [xx(7), xx(1)*xx(8), xx(2)*xx(9); xx(3)*xx(7), xx(8), xx(4)*xx(9); xx(5)*xx(7), xx(6)*xx(8),xx(9)]*[gyro_x - 0.008968; gyro_y + 0.010061; gyro_z + 0.015494];
    
    gyM = [xx(1), xx(2), xx(3); xx(4), xx(5), xx(6); xx(7), xx(8),xx(9)]*[gyro_x - 0.008968; gyro_y + 0.010061; gyro_z + 0.015494];
    
    gyroX = gyM(1);
    gyroY = gyM(2);
    gyroZ = gyM(3);
end