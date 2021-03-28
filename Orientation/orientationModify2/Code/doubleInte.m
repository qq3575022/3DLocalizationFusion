function pos = doubleInte(t, a)

vx  = NaN(1,length(a));  vx(1) = 0;
pos = NaN(1,length(a)); pos(1) = 0;

for i = 2:1:length(a)
    
    vx(i)  = vx(i-1) + a(i-1)*(t(i)-t(i-1));
    pos(i) = vx(i)*(t(i)-t(i-1)) + 0.5*a(i)*(t(i)-t(i-1))^2 + pos(i-1);
    
end


end