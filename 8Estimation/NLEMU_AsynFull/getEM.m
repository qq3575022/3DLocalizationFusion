function E = getEM(y,x)

HN = NaN(3,1);

HN(1:3)  = getHxkM(x);

E = y - HN;

end