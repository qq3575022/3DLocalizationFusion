function y = getyNP(r1,r2,r3,r4,k,N)

y = NaN(4*N,1);

for i = 1:1:N
y(1+4*(i-1):4*i) = [r1(k+i-1);r2(k+i-1);r3(k+i-1);r4(k+i-1)];
end

end