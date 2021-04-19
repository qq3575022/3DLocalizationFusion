function F = getF(time, m, m2, accT, i)
 
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
    
F = I +Add*(time(m2) - time(m));

end