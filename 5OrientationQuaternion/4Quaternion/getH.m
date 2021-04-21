function H = getH(x, index)

if index == 1
    H = getHxkPA(x); 
    
elseif index == 2
    H = getHxkA(x); 
elseif index == 3
    H = getHxkP(x); 
    
end

end