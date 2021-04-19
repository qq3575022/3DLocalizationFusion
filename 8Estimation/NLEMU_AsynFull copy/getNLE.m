function E = getNLE(y, x, N, index, len, time, m)

HN = NaN(length(y),1);
Gain = eye(length(x));

for i = 1:1:N

     Gain = getF(time, m, m+i-1);

    
    if index(i) == 1
        HN(len(i)+1:len(i+1)) = getHxkMGA(Gain*x);

    elseif index(i) == 2
        HN(len(i)+1:len(i+1)) = getHxkGA(Gain*x);

    elseif index(i) == 3
        HN(len(i)+1:len(i+1)) = getHxkMG(Gain*x);

    elseif index(i) == 4
        HN(len(i)+1:len(i+1)) = getHxkMA(Gain*x);

    elseif index(i) == 5
        HN(len(i)+1:len(i+1)) = getHxkA(Gain*x);

    elseif index(i) == 6
        HN(len(i)+1:len(i+1)) = getHxkG(Gain*x);

    elseif index(i) == 7
        HN(len(i)+1:len(i+1)) = getHxkM(Gain*x);
    end

end

E = (y - HN).^2;

end