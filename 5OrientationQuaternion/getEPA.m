function E = getEPA(y,x)

HN = NaN(6,1);

HN(1:6)  = getHxkPA(x);

E = y - HN;

end