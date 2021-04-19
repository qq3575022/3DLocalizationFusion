function E = getEA(y,x,a)

HN = NaN(3,1);

HN(1:3)  = getHxkA(x,a);

E = y - HN;

end