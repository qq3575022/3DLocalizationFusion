function [phi,phi_gt,AXY,XX,Tmag,tdmag,Tacc,tdacc,T3,td3] = groundtruth3D()
h = 1e-5;   tc = 0:h:2.1-h;   

w_a = zeros(1,length(tc));    w_a(1) = -0.3*pi;   w_a(length(w_a)) = 0.3*pi;
w_v = pi*ones(1,length(tc));   
w   = zeros(1,length(tc));

for n = 2:length(tc)
    w(n) = w(n-1) +w_v(n)*h;
end

% change time stamp of arrival
tdmag = zeros(1,1409);        tdmag(2:705) = 0.00138:0.00138:0.97152;   tdmag(706:1409) = 0.97298:0.00146:1.99936; % Time Stamp
Tmag  = zeros(1,1409);         Tmag(1:705) = 0.00138;                    Tmag(706:1409) = 0.00146;                 % Sampling Period T

tdacc = 0:0.00142:1.99936;
Tacc  = 0.00142*ones(1,1409); 

T3  = 1.0* ones(1,2818); 
td3 = 1.0*zeros(1,2818);

o = 1;
p = 1;

for iii = 2:2818
    
    if o <= 1409 && p <= 1409 && tdmag(o) < tdacc(p)
        td3(iii) = tdmag(o);
        T3(iii) = td3(iii) - td3(iii - 1);
        o = o + 1;
    else
        if p <= 1409
        td3(iii) = tdacc(p);
        T3(iii) = td3(iii) - td3(iii - 1);
        p = p + 1;
        end
    end
    
end


%-----------------------------------------------------------------------------------------------------------------------------------------------
% 1. Sampling from loaded groudtruth in continuous time - Gyro Mag Simulation
W = NaN(1,length(tdmag)); W_V = NaN(1,length(tdmag)); W_A = NaN(1,length(tdmag));  
W_A(1) = w_a(1);     W_V(1) = w_v(1);     W(1)   = w(1);      

k = 1;
for n = 137:97155
   if mod(n-138,138) == 0
        k = k + 1;
        W_A(k) = w_a(n);                                       % Sampling angular acceleration
        W_V(k) = w_v(n);                                       % Sampling angular velocity
        W(1,k) = w(n);                                         % Sampling angle
    end   
end

for n = 97158:201000
   if mod(n-97298,146) == 0
        k = k + 1;
        W_A(k) = w_a(n);                                      % Sampling angular acceleration
        W_V(k) = w_v(n);                                      % Sampling angular velocity
        W(1,k) = w(n);                                        % Sampling angle
    end   
end

% 1.2 Sampling from loaded groudtruth in continuous time - Gyro Mag Simulation      
W_gt = NaN(1,2*length(tdmag)); W_V_gt = NaN(1,2*length(tdmag)); 
W_gt(1)   = w(1);   W_V_gt(1) = w_v(1);       

k = 1;
for n = 137:97155
   if mod(n-138,138) == 0 || mod(n-138,142) == 0
        k = k + 1;
        W_V_gt(k) = w_v(n);                                       % Sampling angular velocity
        W_gt(1,k) = w(n);                                         % Sampling angle
    end   
end

for n = 97158:201000%190991
   if mod(n-97298,142) == 0 || mod(n-97298,146) == 0 
        k = k + 1;
        W_V_gt(k) = w_v(n);                                      % Sampling angular velocity
        W_gt(1,k) = w(n);                                        % Sampling angle
    end   
end

phi_gt = NaN(6,length(W_V_gt)); 
phi_gt(1,:) = zeros(1,length(W_V_gt)); phi_gt(2,:) = zeros(1,length(W_V_gt)); 
phi_gt(3,:) = zeros(1,length(W_V_gt)); phi_gt(4,:) = zeros(1,length(W_V_gt)); 
phi_gt(5,:) = W_gt;                    phi_gt(6,:) = W_V_gt;

%------------------------------------------------------------------------------------------------------------------------------------------------
% 2. Get position x,y,z; Velocity xdot, ydot, zdot; Acceleration xdotdot, ydotdot, zdotdot - Accelerometer Simu
x = NaN(3,length(tc)); y = NaN(3,length(tc)); z = NaN(3,length(tc)); XX = NaN(9,length(td3)); AXY = NaN(3,length(tdmag)); 

x(1,1) = -0.3*sin(w(1)) + 1.95; x(2,1) = 0; x(3,1) = 0;
y(1,1) =  0.3*cos(w(1)) + 1.51; y(2,1) = 0; y(3,1) = 0;
z(1,1) =  0;                    z(2,1) = 0; z(3,1) = 0;

k = 1;k2 = 1;
for n = 2:201390
    x(1,n) = -0.3*sin(w(n)) + 1.95;                          % Continuous time position - x
    x(2,n) = -0.3*cos(w(n))*w_v(n);                          % Continuous time velocity - xdot
    x(3,n) =  0.3*sin(w(n))*w_v(n)^2 - 0.3*cos(w(n))*w_a(n); % Continuous time acceleration - xdotdot
    
    y(1,n) =  0.3*cos(w(n)) + 1.51;                          % Continuous time position - y
    y(2,n) = -0.3*sin(w(n))*w_v(n);                          % Continuous time velocity - ydot
    y(3,n) = -0.3*cos(w(n))*w_v(n)^2 - 0.3*sin(w(n))*w_a(n); % Continuous time acceleration - ydotdot
    
    z(1,n) = 0;                                              % Continuous time position - z
    z(2,n) = 0;                                              % Continuous time velocity - zdot
    z(3,n) = 0;                                              % Continuous time acceleration - zdotdot
    
    H = getHxkPVA([x(1,n),x(2,n),x(3,n),y(1,n),y(2,n),y(3,n),z(1,n),z(2,n),z(3,n),0,0,0,0,w(n),w_v(n)]);
    
    if (137<n)&&(n<97155)
        if mod(n-138,138) == 0 || mod(n-138,142) == 0
            k = k+1;
            
            XX(1,k) = x(1,n);                                % sampled x
            XX(2,k) = x(2,n);                                % sampled xdot
            XX(3,k) = x(3,n);                                % sampled xdotdot

            XX(4,k) = y(1,n);                                % sampled y
            XX(5,k) = y(2,n);                                % sampled ydot
            XX(6,k) = y(3,n);                                % sampled ydotdot

            XX(7,k) = z(1,n);                                % sampled z
            XX(8,k) = z(2,n);                                % sampled zdot
            XX(9,k) = z(3,n);                                % sampled zdotdot

        end
     else 
        if mod(n-97298,142) == 0 || mod(n-97298,146) == 0 
            k = k+1;      
            XX(1,k) = x(1,n); 
            XX(2,k) = x(2,n); 
            XX(3,k) = x(3,n);

            XX(4,k) = y(1,n); 
            XX(5,k) = y(2,n); 
            XX(6,k) = y(3,n);

            XX(7,k) = z(1,n);
            XX(8,k) = z(2,n);
            XX(9,k) = z(3,n);

        end
    end
    
    if mod(n,142) == 0
        k2 = k2+1;      
        AXY(1,k2) = H(1);
        AXY(2,k2) = H(2);
        AXY(3,k2) = H(3);
    end
    
end
    
XX(1,1)  = x(1,142);  XX(2,1)  = x(2,142);  XX(3,1)  = x(3,142);
XX(4,1)  = y(1,142);  XX(5,1)  = y(2,142);  XX(6,1)  = y(3,142);
XX(7,1)  = z(1,142);  XX(8,1)  = z(2,142);  XX(9,1)  = z(3,142);
AXY(1,1) = 0;         AXY(2,1) = 0;         AXY(3,1) = 0;


%------------------------------------------------------------------------------------------------------------------------------------------------
% 3. Generate Simulated Input Sensor Data
% magtometer and gyroscope measurements roll pitch yaw: [\phi \phi_dot \theta \theta_dot \psi \psi_dot]

phi = NaN(6,length(W_V)); 
phi(1,:) = zeros(1,length(W_V)); phi(2,:) = zeros(1,length(W_V)); 
phi(3,:) = zeros(1,length(W_V)); phi(4,:) = zeros(1,length(W_V)); 
phi(5,:) = W;                    phi(6,:) = W_V;

% Add noise
phi = awgn(phi,20);
AXY = awgn(AXY,20);

end