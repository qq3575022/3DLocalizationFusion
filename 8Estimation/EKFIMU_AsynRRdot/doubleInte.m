function pos = doubleInte(T, a)

vx  = NaN(1,length(a));  vx(1) = 0;
pos = NaN(1,length(a)); pos(1) = 0;

for i = 2:1:length(a)
    vx(i)  = vx(i-1) + a(i-1)*T(i);
    pos(i) = vx(i)*T(i) + 0.5*a(i)*T(i)^2 + pos(i-1);

end


end