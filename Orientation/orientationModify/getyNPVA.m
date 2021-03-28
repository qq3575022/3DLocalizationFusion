function [y, index, n1] = getyNPVA(tdgyro, tdacc, gyro, AXY, n1)
% 
% td3m = tdgyro
% tdaccn1 = tdacc(n1)

thres = 100;

if (n1 <= length(AXY)) & (abs(round(tdgyro*1e5)-round(tdacc(n1)*1e5)) < thres) 
    
    y = NaN(3,1);
    y(1:3) = [AXY(1,n1); AXY(2,n1); AXY(3,n1)];
    
    n1 = n1 + 1;
    
    index = 1;
      

else
    y = NaN(3,1);index = -1;
end 


end