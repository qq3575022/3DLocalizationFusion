function E = getNLE(y, x, N, time, m, accT, i)

HN = NaN(length(y),1);
Gain = eye(length(x));

for i = 1:1:N

     Gain = getF(time, m, m+i-1, accT, i);
     HN(1+9*(i-1):9*i) = getHxkMGA(Gain*x);

end

E = (y - HN).^2;

end