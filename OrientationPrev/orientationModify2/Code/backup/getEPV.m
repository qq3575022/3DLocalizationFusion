function E = getEPV(y,x)

HN = NaN(6,1);

HN(1:6)  = getHxkPV(x);

E = y - HN;

end