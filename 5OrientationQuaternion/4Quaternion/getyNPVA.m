function [y, index, n1, n2] = getyNPVA(tdgyro, tdacc, tdmag, mag, gyro, AXY, n1, n2)

% td3m = td3m
% tdaccn1 = tdacc(n1)
% tdgyron2 = tdmag(n2)

thres = 10;

if n1 <= length(AXY) & abs(round(tdgyro*1e5)-round(tdacc(n1)*1e5)) < thres & n2 <= length(mag) & abs(round(tdgyro*1e5) - round(tdmag(n2)*1e5)) < thres
    
    y = NaN(6,1);
    y(1:6) = [AXY(1,n1); AXY(2,n1); AXY(3,n1); mag(1,n2); mag(2,n2); mag(3,n2)];
    
    n1 = n1 + 1;
    n2 = n2 + 1;
    
    index = 1;
        
elseif n1 <= length(AXY) & abs(round(tdgyro*1e5)-round(tdacc(n1)*1e5)) < thres
    
    y = NaN(3,1);
    y(1:3) = [AXY(1,n1); AXY(2,n1); AXY(3,n1)];
    
    n1 = n1 + 1;
    
    index = 2;
        
elseif n2 <= length(mag) & abs(round(tdgyro*1e5) - round(tdmag(n2)*1e5)) < thres
    
    y = NaN(3,1);
    y(1:3) = [mag(1,n3); mag(2,n3); mag(3,n3)];
       
    n2 = n2 + 1;
    
    index = 3;
        

else
    y = NaN(3,1);index = -1;
end 


end