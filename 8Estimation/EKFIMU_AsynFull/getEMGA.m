function E = getEMGA(y,x)

HN = NaN(9,1);

HN(1:9)  = getHxkMGA(x);

E = y - HN;

end