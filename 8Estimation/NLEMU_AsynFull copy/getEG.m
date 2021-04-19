function E = getEG(y,x)

HN = NaN(3,1);

HN(1:3)  = getHxkG(x);

E = y - HN;

end