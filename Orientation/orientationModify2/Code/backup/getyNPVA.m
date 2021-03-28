function [y, index, n1, n2, n3] = getyNPVA(td3m, tdacc, tdgyro, tdmag, mag, gyro, AXY, n1, n2, n3)

% td3m = td3m
% tdaccn1 = tdacc(n1)
% tdgyron2 = tdgyro(n2)
% tymagn3 = tdmag(n3)

if n1 <= length(AXY) & abs(round(td3m*1e5)-round(tdacc(n1)*1e5)) <= 10 & n2 <= length(gyro) & abs(round(td3m*1e5) - round(tdgyro(n2)*1e5)) <= 10 & n3 <= length(mag) & abs(round(td3m*1e5) - round(tdmag(n3)*1e5)) <= 10
    
    y = NaN(9,1);
    y(1:9) = [AXY(1,n1); AXY(2,n1); AXY(3,n1); mag(1,n3); gyro(1,n2); mag(2,n3); gyro(2,n2); mag(3,n3); gyro(3,n2)];
    
    n1 = n1 + 1;
    n2 = n2 + 1;
    n3 = n3 + 1;
    
    index = 1;
        
elseif n1 <= length(AXY) & abs(round(td3m*1e5)-round(tdacc(n1)*1e5)) <= 10 & n2 <= length(gyro) &  abs(round(td3m*1e5) - round(tdgyro(n2)*1e5)) <= 10
    
    y = NaN(6,1);
    y(1:6) = [AXY(1,n1); AXY(2,n1); AXY(3,n1); gyro(1,n2); gyro(2,n2); gyro(3,n2)];
    
    n1 = n1 + 1;
    n2 = n2 + 1;
    
    index = 2;
        
elseif n1 <= length(AXY) & abs(round(td3m*1e5)-round(tdacc(n1)*1e5)) <= 10 & n3 <= length(mag) & abs(round(td3m*1e5) - round(tdmag(n3)*1e5)) <= 10
    
    y = NaN(6,1);
    y(1:6) = [AXY(1,n1); AXY(2,n1); AXY(3,n1); mag(1,n3); mag(2,n3); mag(3,n3)];
       
    n1 = n1 + 1;
    n3 = n3 + 1;
    
    index = 3;
        
elseif n2 <= length(gyro) & abs(round(td3m*1e5) - round(tdgyro(n2)*1e5)) <= 10 & n3 <= length(mag) & abs(round(td3m*1e5) - round(tdmag(n3)*1e5)) <= 10
    
    y = NaN(6,1);
    y(1:6) = [mag(1,n3); gyro(1,n2); mag(2,n3); gyro(2,n2); mag(3,n3); gyro(3,n2)];
    
    n2 = n2 + 1;
    n3 = n3 + 1;
    
    index = 4;
    
elseif n1 <= length(AXY) & abs(round(td3m*1e5)-round(tdacc(n1)*1e5)) <= 10
    
    y = NaN(3,1);
    y(1:3) = [AXY(1,n1); AXY(2,n1); AXY(3,n1);];
    
    n1 = n1 + 1;
    
    index = 5;

elseif n2 <= length(gyro) & abs(round(td3m*1e5) - round(tdgyro(n2)*1e5)) <= 10
    y = NaN(3,1);
    y(1:3) = [gyro(1,n2); gyro(2,n2); gyro(3,n2)];
           
    n2 = n2 + 1;
    
    index = 6;
        
elseif n3 <= length(mag) & abs(round(td3m*1e5) - round(tdmag(n3)*1e5)) <= 10
    y = NaN(3,1);
    y(1:3) = [mag(1,n3); mag(2,n3); mag(3,n3)];
        
    n3 = n3 + 1;
    
    index = 7;
else
    y = NaN(3,1);index = -1;
end 


end