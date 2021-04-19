function E = getEP(y,x,g)

E = 0;

for i = 1:1:length(y)
      
   E = E + ((x(4)*(y(i,1) + x(7)) + x(1)*x(5)*(y(i,2) + x(8)) + x(2)*x(6)*(y(i,3) + x(9)))^2 + (x(5)*(y(i,2) + x(8)) + x(3)*x(6)*(y(i,3) + x(9)))^2 + (x(6)*(y(i,3) + x(9)))^2  - g^2)^2 ;

end


end