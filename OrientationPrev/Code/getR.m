function R = getR(RR, index)


if index == 1
    R = diag([RR(1), RR(2), RR(3), RR(4), RR(5), RR(6), RR(7), RR(8), RR(9)]); 

% 2-4
elseif index == 2 
    R = diag([RR(1), RR(2), RR(3), RR(5), RR(7), RR(9)]); 
elseif index == 3 
    R = diag([RR(1), RR(2), RR(3), RR(4), RR(6), RR(8)]); 
elseif index == 4
    R = diag([RR(4), RR(5), RR(6), RR(7), RR(8), RR(9)]); 

% 5-7    
elseif index == 5
    R = diag([RR(1), RR(2), RR(3)]); 
elseif index == 6
    R = diag([RR(5), RR(7), RR(9)]); 
elseif index == 7
    R = diag([RR(4), RR(6), RR(8)]); 
else
    R = nan;
end

end