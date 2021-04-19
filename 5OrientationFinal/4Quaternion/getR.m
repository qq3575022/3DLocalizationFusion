function R = getR(RR, index)


if index == 1
    R = diag([RR(1), RR(2), RR(3), RR(4), RR(5), RR(6)]); 

% 2-4
elseif index == 2 
    R = diag([RR(1), RR(2), RR(3)]); 
elseif index == 3 
    R = diag([RR(4), RR(5), RR(6)]); 
else
    R = nan;
end

end