function E = getEMG(y,x)

HN = NaN(6,1);

HN(1:6)  = getHxkMG(x);

E = y - HN;

end