function H = getH(x, index)

if index == 1
    H = getHxkMGA(x);

elseif index == 2
    H = getHxkGA(x);

elseif index == 3
    H = getHxkMG(x); 
    
elseif index == 4
    H = getHxkMA(x);
    
elseif index == 5
    H = getHxkA(x);

elseif index == 6
    H = getHxkG(x); 
    
elseif index == 7    
    H = getHxkM(x); 
end

end