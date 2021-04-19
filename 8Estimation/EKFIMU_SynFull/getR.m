function R = getR(RR, index)
%[mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z];
%var(accxZ), var(accyZ), var(acczZ), var(magxZ), var(gyroxZ), var(magyZ), var(gyroyZ), var(magzZ), var(gyrozZ)]


%         E = getEA(y,x-dx*e);
%         E = getEG(y,x-dx*e);
%         E = getEM(y,x-dx*e);


if index == 1
    R = diag([RR(1), RR(2), RR(3), RR(4), RR(5), RR(6), RR(7), RR(8), RR(9)]); 

elseif index == 2
    R = diag([RR(1), RR(2), RR(3), RR(5), RR(7), RR(9)]); 

elseif index == 3
    R = diag([RR(4), RR(5), RR(6), RR(7), RR(8), RR(9)]); 
    
elseif index == 4
    R = diag([RR(1), RR(2), RR(3), RR(4), RR(6), RR(8)]); 
    
elseif index == 5
    R = diag([RR(1), RR(2), RR(3)]); 

elseif index == 6
    R = diag([RR(5), RR(7), RR(9)]); 
    
elseif index == 7
    R = diag([RR(4), RR(6), RR(8)]); 
else
    R = nan;


   
        
end