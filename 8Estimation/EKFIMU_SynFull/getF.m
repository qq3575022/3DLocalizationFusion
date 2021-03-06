function F = getF(time, m, accT, gyroT, magT, i, j, k)
 
I   =  [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
   
Add  = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    %accT, gyroT, magT, i, j, k
F = I +Add*[0; time(m)-time(m-1); time(m)-time(m-1); 0; time(m)-time(m-1); time(m)-time(m-1); 0; time(m)-time(m-1); time(m)-time(m-1); 0; magT(k-1)-time(m-1); 0; magT(k-1)-time(m-1); 0; magT(k-1)-time(m-1);];

end