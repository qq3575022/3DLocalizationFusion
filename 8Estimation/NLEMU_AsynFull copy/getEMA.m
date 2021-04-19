function E = getEMA(y,x)

HN = NaN(6,1);

HN(1:6)  = getHxkMA(x);

E = y - HN;

end