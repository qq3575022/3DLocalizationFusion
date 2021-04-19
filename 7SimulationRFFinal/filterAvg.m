function [time2, magD3, magD4, magD5, magD6, magD7, phaseD2] = filterAvg(magD1, phaseD1)



%d2=unwrap(phaseD1);
magD2 = zeros(1, round(length(magD1)/100));
magD3 = zeros(1, round(length(magD1)/100));
magD4 = zeros(1, round(length(magD1)/100));
magD5 = zeros(1, round(length(magD1)/100));
magD6 = zeros(1, round(length(magD1)/100));
magD7 = zeros(1, round(length(magD1)/100));

for i = 2:1:length(magD1)/100
     magD2(i) = mean(magD1(100*(i-1)+1:100*i));
    
end

time2 = 0:length(magD2)/(length(magD2)*60):length(magD2)/60 - length(magD2)/(length(magD2)*60);

for i = 1:1:length(magD2) - 100
     magD3(i) = mean(magD2(i:100 + i));
end


for i = 1:1:length(magD3) - 100
     magD4(i) = mean(magD3(i:100 + i));
end

for i = 1:1:length(magD4) - 100
     magD5(i) = mean(magD4(i:100 + i));
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

for i = 1:1:length(magD5) - 100
     magD6(i) = mean(magD5(i:100 + i));
end

for i = 1:1:length(magD6) - 100
     magD7(i) = mean(magD6(i:100 + i));
end

phaseD2 = zeros(1, round(length(phaseD1)/100));

for i = 1:1:length(phaseD1)/100
    phaseD2(i) = mean(phaseD1(100*(i-1)+1:100*i));
end

d3 = diff(phaseD2)/10e-3;

% for i = 2:length(d3)
%     if (d3(i) - d3(i-1)) > 8
%         d3(i) = d3(i-1);
%     end
% end

d4  = zeros(1, length(d3));
d4T = zeros(1, length(d3));

for j = 1:1:length(d4)-6
    %phaseD2(i) = mean(phaseD1(1000*(i-1)+1:1000*i));
    d4(j) = mean(d3(j:5+j));%/10^(-6);
    d4T(j) = j/100;
    
end

end