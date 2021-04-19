function E = getEP(y,x)

E = 0;

for i = 1:1:length(y)
      
   E = E + (y(1,i) - x(1))^2 + (y(2,i) - x(2))^2 + (y(3,i) - x(3))^2 - 10^4;
end


end