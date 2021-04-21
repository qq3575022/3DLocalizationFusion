%function [x, y, z, time] = getPath()
% Feb 14th
time = 0: 0.0001:20;
            
% % get angular velocity w
w_a = zeros(1,length(time));w_v = zeros(1,length(time));w = zeros(1,length(time)); T = mean(diff(time));

w_a(time > 2.8 & time <= 8.6) = 0.125;
w_a(time > 8.6 & time <= 11.5) = 0;
w_a(time > 11.5 & time <= 17.3) = -0.125;

%load('X.mat');load('A.mat');load('AxB.mat');load('AyB.mat')
%w = X(1,1:190991);     w_v = X(2,1:190991);     w_a = A(1,1:190991);  

W = NaN(1,round(length(time)/188));   W_V = NaN(1,round(length(time)/25)); %W_A = NaN(1,length(td));  
t_m = NaN(1,round(length(time)/188)); t_v = NaN(1,round(length(time)/25)); 

km = 1;
kg = 1;

for n = 2:length(w)
    
    w_v(n) = w_v(n-1) + w_a(n-1)*T;
    w(n) = w(n-1) + w_v(n-1)*T + 0.5*w_a(n-1)*T^2;
    
    if mod(n,25) == 0
        t_v(kg) = time(n);
        W_V(kg) = w_v(n);
        kg = kg + 1;

    end   
    
    if mod(n,188) == 0
        t_m(km) = time(n);
        W(km) = w(n);
        km = km + 1;

    end  
    
end

W(1) = w(1); W_V(1) = w_v(1);

figure
subplot(211), plot(t_v, W_V)
subplot(212), plot(t_m, W)

% get position x,y
x = NaN(3,length(w));  y = NaN(3,length(w)); 
ax = NaN(1,length(w));ay = NaN(1,length(w));

XX = NaN(6,round(length(time)/188)+round(length(time)/25)-42); AXY = NaN(2,round(length(time)/25));
t_a = NaN(1,round(length(time)/25)); 
t_x = NaN(1,round(length(time)/188)+round(length(time)/25)-42);

x(1,1) =  -0.3*sin(w(1)) + 0.3;  x(2,1) = 0; x(3,1) = 0;
y(1,1) =   0.3*cos(w(1)) + 0.3;  y(2,1) = 0; y(3,1) = 0;


ka = 1;
kx = 1;

for n = 2:length(w)-1
    x(1,n) =  0.3*sin(-w(n)) + 0.3;
    x(2,n) = -0.3*cos(-w(n))*w_v(n);
    x(3,n) = -0.3*sin(-w(n))*w_v(n)^2 - 0.3*cos(-w(n))*w_a(n);
    
    y(1,n) =  0.3*cos(-w(n)) + 0.3;
    y(2,n) =  0.3*sin(-w(n))*w_v(n);
    y(3,n) = -0.3*cos(-w(n))*w_v(n)^2 + 0.3*sin(-w(n))*w_a(n);
    
    ax(1,n) =  x(3,n)*sin(-w(n))+y(3,n)*cos(-w(n));
    ay(1,n) = -x(3,n)*cos(-w(n))+y(3,n)*sin(-w(n));
    
    if mod(n,25) == 0
       
        XX(1,kx) = x(1,n); 
        XX(2,kx) = x(2,n); 
        XX(3,kx) = x(3,n);
        
        AXY(1,ka) = ax(1,n);
        AXY(2,ka) = ay(1,n);
        
        t_a(ka) = time(n);
        
        XX(4,kx) = y(1,n); 
        XX(5,kx) = y(2,n); 
        XX(6,kx) = y(3,n); 
        
        t_x(kx) = time(n);
        
        ka = ka+1;  
        kx = kx+1;  
        
    end
    
    if mod(n,188) == 0
       
        XX(1,kx) = x(1,n); 
        XX(2,kx) = x(2,n); 
        XX(3,kx) = x(3,n);
        
        
        XX(4,kx) = y(1,n); 
        XX(5,kx) = y(2,n); 
        XX(6,kx) = y(3,n); 
        
        t_x(kx) = time(n);
          
        kx = kx+1;  
        
    end
    
    
end

XX(1,1) = x(1,1); XX(2,1) = 0; XX(3,1) = 0;
XX(4,1) = y(1,1); XX(5,1) = 0; XX(6,1) = 0; 

        
figure
subplot(211), plot(t_a, AXY(1,:))
subplot(212), plot(t_a, AXY(2,:))   

figure
plot(XX(1,:), XX(4,:));
%end