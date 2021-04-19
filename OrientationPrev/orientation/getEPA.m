function E = getEPA(y,x)

HN = NaN(3,1);

HN(1:3)  = getHxkPA(x);

E = y - HN;

end