function E = getEA(y,x)

HN = NaN(3,1);

HN(1:3)  = getHxkA(x);

E = y - HN;

end