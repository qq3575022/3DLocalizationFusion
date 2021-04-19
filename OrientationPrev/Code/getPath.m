function [t_v, T_v, W_V, t_m, T_m, W, t_a, T_a, AXY, t_x, XX] = getPath()
%[phi,phi_gt,AXY,XX,Tmag,tdmag,Tacc,tdacc,T3,td3] 
h = 0.0001;
time = 0: h :20;
            
% % get angular velocity w
w_a = zeros(1,length(time));w_v = zeros(1,length(time));w = zeros(1,length(time)); T = mean(diff(time));

% w_a = zeros(1,length(time));    w_a(1) = -0.03*pi;   w_a(length(w_a)) = 0.03*pi;
% w_v = pi/10*ones(1,length(time));   
% w   = zeros(1,length(time));

% for n = 2:length(time)
%     w(n) = w(n-1) +w_v(n)*h;
% end
%

%%
% w_a(time > 2.8 & time <= 8.6) = 0.125;
% w_a(time > 8.6 & time <= 11.5) = 0;
% w_a(time > 11.5 & time <= 17.3) = -0.125;

W   = NaN(3,round(length(time)/188));   W_V = NaN(3,round(length(time)/25)); %W_A = NaN(1,length(td));  
t_m = NaN(1,round(length(time)/188));   t_v = NaN(1,round(length(time)/25)); 
T_m = NaN(1,round(length(time)/188));   T_v = NaN(1,round(length(time)/25)); 

km = 1;
kg = 1;

for n = 2:length(w)
    
    
    if mod(n,25) == 0
        t_v(kg) = time(n);
        T_v(kg) = 25;
        W_V(1, kg) = 0;
        W_V(2, kg) = 0;
        W_V(3, kg) = w_v(n);
        
        kg = kg + 1;

    end   
    
    if mod(n,188) == 0
        t_m(km) = time(n);
        T_m(km) = 25;
        W(1, km) = 0;
        W(2, km) = 0;
        W(3, km) = w(n);

        km = km + 1;

    end  
    
end

W(1) = w(1); W_V(1) = w_v(1);

% get position x,y
x  = NaN(3,length(w));  y = NaN(3,length(w)); 
ax = NaN(1,length(w)); ay = NaN(1,length(w));

XX  = NaN(6,round(length(time)/188)+round(length(time)/25)-42); AXY = NaN(3,round(length(time)/25));
t_a = NaN(1,round(length(time)/25)); 
T_a = NaN(1,round(length(time)/25)); 
t_x = NaN(1,round(length(time)/188)+round(length(time)/25)-42);
T_x = NaN(1,round(length(time)/188)+round(length(time)/25)-43);

x(1,1) =  -0.3*sin(w(1)) + 0.3;  x(2,1) = 0; x(3,1) = 0;
y(1,1) =   0.3*cos(w(1)) - 0.3;  y(2,1) = 0; y(3,1) = 0;


ka = 1;
kx = 1;

for n = 2:length(w)
    x(1,n) = -0.3*sin(w(n));% + 1.95;                          % Continuous time position - x
    x(2,n) = -0.3*cos(w(n))*w_v(n) + 0.3*cos(w(1))*w_v(1);                          % Continuous time velocity - xdot
    x(3,n) =  0.3*sin(w(n))*w_v(n)^2 - 0.3*cos(w(n))*w_a(n); % Continuous time acceleration - xdotdot
    
    y(1,n) =  0.3*cos(w(n)) - 0.3;% + 1.51;                          % Continuous time position - y
    y(2,n) = -0.3*sin(w(n))*w_v(n);                          % Continuous time velocity - ydot
    y(3,n) = -0.3*cos(w(n))*w_v(n)^2 - 0.3*sin(w(n))*w_a(n); % Continuous time acceleration - ydotdot
    
    
    ax(1,n) =  x(3,n)*cos(w(n))-y(3,n)*sin(w(n));            % Continuous time body axis acceleration - ax
    ay(1,n) =  x(3,n)*sin(w(n))+y(3,n)*cos(w(n));            % Continuous time body axis acceleration - ay

%     
    if mod(n,25) == 0
       
        XX(1,kx) = x(1,n); 
        XX(2,kx) = x(2,n); 
        XX(3,kx) = x(3,n);
        
        AXY(1,ka) = ax(1,n);
        AXY(2,ka) = ay(1,n);
        AXY(3,ka) = 0;
        
        t_a(ka) = time(n);
        T_a(ka) = 25;
        
        XX(4,kx) = y(1,n); 
        XX(5,kx) = y(2,n); 
        XX(6,kx) = y(3,n); 
        
        t_x(kx) = time(n);
        
        if kx == 1
            T_x(kx) = 0;
        else
            T_x(kx) = time(n) - t_x(kx-1);
        end
        
        
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
        
        if kx == 1
            T_x(kx) = 0;
        else
            T_x(kx) = time(n) - t_x(kx-1);
        end
        
        kx = kx+1;  
        
    end
    
    
end

% 6 data
% AXY = zeros(3, length(t_a));
% AXY(:,1:round(length(t_a)/6)) = [0;0;9.793]*ones(1,round(length(t_a)/6));
% AXY(:,round(length(t_a)/6)+1:2*round(length(t_a)/6)) = [0;9.793;0]*ones(1,round(length(t_a)/6));
% AXY(:,2*round(length(t_a)/6)+1:3*round(length(t_a)/6)) = [9.793;0;0]*ones(1,round(length(t_a)/6));
% 
% AXY(:,3*round(length(t_a)/6)+1:4*round(length(t_a)/6)) = [9.793*0.5;0*0.5*sqrt(3);0]*ones(1,round(length(t_a)/6));
% AXY(:,4*round(length(t_a)/6)+1:5*round(length(t_a)/6)) = [9.793*0.5;0;0*0.5*sqrt(3)]*ones(1,round(length(t_a)/6));
% AXY(:,5*round(length(t_a)/6)+1:6*round(length(t_a)/6)) = [0;9.793*0.5;0*0.5*sqrt(3)]*ones(1,round(length(t_a)/6));

% 12 data
% AXY = zeros(3, length(t_a));
% AXY(:,1:round(length(t_a)/12)) = [0;0;9.793]*ones(1,round(length(t_a)/12));
% AXY(:,round(length(t_a)/12)+1:2*round(length(t_a)/12))   = [0;9.793;0]*ones(1,round(length(t_a)/12));
% AXY(:,2*round(length(t_a)/12)+1:3*round(length(t_a)/12)) = [9.793;0;0]*ones(1,round(length(t_a)/12));
% 
% AXY(:,3*round(length(t_a)/12)+1:4*round(length(t_a)/12)) = [9.793*0.5;9.793*0.5*sqrt(3);0]*ones(1,round(length(t_a)/12));
% AXY(:,4*round(length(t_a)/12)+1:5*round(length(t_a)/12)) = [9.793*0.5;0;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/12));
% AXY(:,5*round(length(t_a)/12)+1:6*round(length(t_a)/12)) = [0;9.793*0.5;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/12));
% 
% AXY(:,6*round(length(t_a)/12)+1:7*round(length(t_a)/12)) = [9.793*0.5*sqrt(2);9.793*0.5*sqrt(2);0]*ones(1,round(length(t_a)/12));
% AXY(:,7*round(length(t_a)/12)+1:8*round(length(t_a)/12)) = [9.793*0.5*sqrt(2);0;9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/12));
% AXY(:,8*round(length(t_a)/12)+1:9*round(length(t_a)/12)) = [0;9.793*0.5*sqrt(2);9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/12));
% 
% AXY(:,9*round(length(t_a)/12)+1:10*round(length(t_a)/12)) = [9.793*0.5*sqrt(3);9.793*0.5;0]*ones(1,round(length(t_a)/12));
% AXY(:,10*round(length(t_a)/12)+1:11*round(length(t_a)/12)) = [9.793*0.5*sqrt(3);0;9.793*0.5]*ones(1,round(length(t_a)/12));
% AXY(:,11*round(length(t_a)/12)+1:length(t_a)) = [0;9.793*0.5*sqrt(3);9.793*0.5]*ones(1,length(t_a) - 11*round(length(t_a)/12));

% n = 18;
% AXY = zeros(3, length(t_a));
% AXY(:,1:round(length(t_a)/n)) = [0;0;9.793]*ones(1,round(length(t_a)/n));
% AXY(:,round(length(t_a)/n)+1:2*round(length(t_a)/n))   = [0;9.793;0]*ones(1,round(length(t_a)/n));
% AXY(:,2*round(length(t_a)/n)+1:3*round(length(t_a)/n)) = [9.793;0;0]*ones(1,round(length(t_a)/n));
% 
% AXY(:,3*round(length(t_a)/n)+1:4*round(length(t_a)/n)) = [9.793*0.5;9.793*0.5*sqrt(3);0]*ones(1,round(length(t_a)/n));
% AXY(:,4*round(length(t_a)/n)+1:5*round(length(t_a)/n)) = [9.793*0.5;0;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/n));
% AXY(:,5*round(length(t_a)/n)+1:6*round(length(t_a)/n)) = [0;9.793*0.5;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/n));
% 
% AXY(:,6*round(length(t_a)/n)+1:7*round(length(t_a)/n)) = [9.793*0.5*sqrt(2);9.793*0.5*sqrt(2);0]*ones(1,round(length(t_a)/n));
% AXY(:,7*round(length(t_a)/n)+1:8*round(length(t_a)/n)) = [9.793*0.5*sqrt(2);0;9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/n));
% AXY(:,8*round(length(t_a)/n)+1:9*round(length(t_a)/n)) = [0;9.793*0.5*sqrt(2);9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/n));
% 
% AXY(:, 9*round(length(t_a)/n)+1:10*round(length(t_a)/n)) = [9.793*0.5*sqrt(3);9.793*0.5;0]*ones(1,round(length(t_a)/n));
% AXY(:,10*round(length(t_a)/n)+1:11*round(length(t_a)/n)) = [9.793*0.5*sqrt(3);0;9.793*0.5]*ones(1,round(length(t_a)/n));
% AXY(:,11*round(length(t_a)/n)+1:12*round(length(t_a)/n)) = [0;9.793*0.5*sqrt(3);9.793*0.5]*ones(1,round(length(t_a)/n));
% 
% AXY(:,12*round(length(t_a)/n)+1:13*round(length(t_a)/n)) = [9.793*sin(5*pi/6);9.793*cos(5*pi/6);0]*ones(1,round(length(t_a)/n));
% AXY(:,13*round(length(t_a)/n)+1:14*round(length(t_a)/n)) = [9.793*sin(5*pi/6);0;9.793*cos(5*pi/6)]*ones(1,round(length(t_a)/n));
% AXY(:,14*round(length(t_a)/n)+1:15*round(length(t_a)/n)) = [0;9.793*sin(5*pi/6);9.793*cos(5*pi/6)]*ones(1,round(length(t_a)/n));
% 
% AXY(:,15*round(length(t_a)/n)+1:16*round(length(t_a)/n)) = [9.793*sin(2*pi/3);9.793*cos(2*pi/3);0]*ones(1,round(length(t_a)/n));
% AXY(:,16*round(length(t_a)/n)+1:17*round(length(t_a)/n)) = [9.793*sin(2*pi/3);0;9.793*cos(2*pi/3)]*ones(1,round(length(t_a)/n));
% 
% AXY(:,17*round(length(t_a)/n)+1:length(t_a)) = [0;9.793*sin(2*pi/3);9.793*cos(2*pi/3)]*ones(1,length(t_a) - 17*round(length(t_a)/n));


n = 24;
AXY = zeros(3, length(t_a));
AXY(:,1:round(length(t_a)/n)) = [0;0;9.793]*ones(1,round(length(t_a)/n));
AXY(:,round(length(t_a)/n)+1:2*round(length(t_a)/n))   = [0;9.793;0]*ones(1,round(length(t_a)/n));
AXY(:,2*round(length(t_a)/n)+1:3*round(length(t_a)/n)) = [9.793;0;0]*ones(1,round(length(t_a)/n));

AXY(:,3*round(length(t_a)/n)+1:4*round(length(t_a)/n)) = [9.793*0.5;9.793*0.5*sqrt(3);0]*ones(1,round(length(t_a)/n));
AXY(:,4*round(length(t_a)/n)+1:5*round(length(t_a)/n)) = [9.793*0.5;0;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/n));
AXY(:,5*round(length(t_a)/n)+1:6*round(length(t_a)/n)) = [0;9.793*0.5;9.793*0.5*sqrt(3)]*ones(1,round(length(t_a)/n));

AXY(:,6*round(length(t_a)/n)+1:7*round(length(t_a)/n)) = [9.793*0.5*sqrt(2);9.793*0.5*sqrt(2);0]*ones(1,round(length(t_a)/n));
AXY(:,7*round(length(t_a)/n)+1:8*round(length(t_a)/n)) = [9.793*0.5*sqrt(2);0;9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/n));
AXY(:,8*round(length(t_a)/n)+1:9*round(length(t_a)/n)) = [0;9.793*0.5*sqrt(2);9.793*0.5*sqrt(2)]*ones(1,round(length(t_a)/n));

AXY(:, 9*round(length(t_a)/n)+1:10*round(length(t_a)/n)) = [9.793*0.5*sqrt(3);9.793*0.5;0]*ones(1,round(length(t_a)/n));
AXY(:,10*round(length(t_a)/n)+1:11*round(length(t_a)/n)) = [9.793*0.5*sqrt(3);0;9.793*0.5]*ones(1,round(length(t_a)/n));
AXY(:,11*round(length(t_a)/n)+1:12*round(length(t_a)/n)) = [0;9.793*0.5*sqrt(3);9.793*0.5]*ones(1,round(length(t_a)/n));

AXY(:,12*round(length(t_a)/n)+1:13*round(length(t_a)/n)) = [9.793*sin(5*pi/6);9.793*cos(5*pi/6);0]*ones(1,round(length(t_a)/n));
AXY(:,13*round(length(t_a)/n)+1:14*round(length(t_a)/n)) = [9.793*sin(5*pi/6);0;9.793*cos(5*pi/6)]*ones(1,round(length(t_a)/n));
AXY(:,14*round(length(t_a)/n)+1:15*round(length(t_a)/n)) = [0;9.793*sin(5*pi/6);9.793*cos(5*pi/6)]*ones(1,round(length(t_a)/n));

AXY(:,15*round(length(t_a)/n)+1:16*round(length(t_a)/n)) = [9.793*sin(2*pi/3);9.793*cos(2*pi/3);0]*ones(1,round(length(t_a)/n));
AXY(:,16*round(length(t_a)/n)+1:17*round(length(t_a)/n)) = [9.793*sin(2*pi/3);0;9.793*cos(2*pi/3)]*ones(1,round(length(t_a)/n));
AXY(:,17*round(length(t_a)/n)+1:18*round(length(t_a)/n)) = [0;9.793*sin(2*pi/3);9.793*cos(2*pi/3)]*ones(1,round(length(t_a)/n));

AXY(:,18*round(length(t_a)/n)+1:19*round(length(t_a)/n)) = [9.793*sin(pi/8);9.793*cos(pi/8);0]*ones(1,round(length(t_a)/n));
AXY(:,19*round(length(t_a)/n)+1:20*round(length(t_a)/n)) = [9.793*sin(pi/8);0;9.793*cos(pi/8)]*ones(1,round(length(t_a)/n));
AXY(:,20*round(length(t_a)/n)+1:21*round(length(t_a)/n)) = [0;9.793*sin(pi/8);9.793*cos(pi/8)]*ones(1,round(length(t_a)/n));

AXY(:,21*round(length(t_a)/n)+1:22*round(length(t_a)/n)) = [9.793*sin(7*pi/8);9.793*cos(7*pi/8);0]*ones(1,round(length(t_a)/n));
AXY(:,22*round(length(t_a)/n)+1:23*round(length(t_a)/n)) = [9.793*sin(7*pi/8);0;9.793*cos(7*pi/8)]*ones(1,round(length(t_a)/n));
AXY(:,23*round(length(t_a)/n)+1:length(t_a)) = [0;9.793*sin(7*pi/8);9.793*cos(7*pi/8)]*ones(1,length(t_a) - 23*round(length(t_a)/n));

%%
stdA = 0.001;
AXY = AXY +stdA*randn(3, length(t_a));


XX(1,:) = doubleInte(t_x, XX(3,:));

XX(1,1) = x(1,1); XX(2,1) = 0; XX(3,1) = 0;
XX(4,1) = y(1,1); XX(5,1) = 0; XX(6,1) = 0;

end