function E = getEGA(y,x)

HN = NaN(6,1);

HN(1:6)  = getHxkGA(x);

E = y - HN;

end