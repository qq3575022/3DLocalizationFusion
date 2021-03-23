function y = getyNPVA(m, AXY)
    
    y = NaN(3,1);
    y(1:3) = [AXY(1,m); AXY(2,m); AXY(3,m)];
    

end