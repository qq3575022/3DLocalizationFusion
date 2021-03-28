function [time2, magD3, magD4, magD5, magD6, magD7] = filterAvg(magD1)



%d2=unwrap(phaseD1);
magD2 = zeros(1, round(length(magD1)/10));
magD3 = zeros(1, round(length(magD1)/10));
magD4 = zeros(1, round(length(magD1)/10));
magD5 = zeros(1, round(length(magD1)/10));
magD6 = zeros(1, round(length(magD1)/10));
magD7 = zeros(1, round(length(magD1)/10));

for i = 2:1:length(magD1)/10
     magD2(i) = mean(magD1(10*(i-1)+1:10*i));
    
end

time2 = 0:length(magD2)/(length(magD2)*60):length(magD2)/60 - length(magD2)/(length(magD2)*60);

for i = 1:1:length(magD2) - 10
     magD3(i) = mean(magD2(i:10 + i));
end


for i = 1:1:length(magD3) - 10
     magD4(i) = mean(magD3(i:10 + i));
end

for i = 1:1:length(magD4) - 10
     magD5(i) = mean(magD4(i:10 + i));
end

% 
for i = 2:1:length(magD5)
%     if magD2(i) > 0.015
%         magD2(i) = magD2(i-1);
%     end
    
    if i > 50000 && magD5(i) > 0.002
        magD5(i) = magD5(i-1);
    end
    if i > 50000 && magD5(i) > 0.002
        magD5(i) = magD5(i-1);
    end
end

for i = 1:1:length(magD5) - 10
     magD6(i) = mean(magD5(i:10 + i));
end

for i = 1:1:length(magD6) - 10
     magD7(i) = mean(magD6(i:10 + i));
end


end