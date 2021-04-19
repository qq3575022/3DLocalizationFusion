function [y, ii, jj, kk] = gety(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m)
    
    [y1, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m);
    ii = i; jj = j; kk = k; 
    [y2, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+1);
    [y3, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+2);
    [y4, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+3);
    [y5, i, j, k] = getyNPVA(mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+4);

    y = [y1;y2;y3;y4;y5];
    
end