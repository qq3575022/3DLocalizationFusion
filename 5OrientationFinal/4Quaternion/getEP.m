function E = getEP(y,x)

HN = NaN(3,1);

HN(1:3)  = getHxkP(x);

E = y - HN;

end