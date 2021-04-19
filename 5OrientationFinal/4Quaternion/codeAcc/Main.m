%%
clc, clear, close all

% Nonlinear - Extended kalman filter
% Load measurement orientation, angular velocity, acc_x, acc_y, acc_z
[mag,gyro,AXY,Tmag,tdmag,Tgyro,tdgyro,Tacc,tdacc,T3,td3] = meas3D();

AXY = [1, 0.0623, 0.0055; 0, 1, -0.0041; 0, 0, 1]*[0.9944, 0, 0; 0, 0.9999, 0; 0, 0, 0.9880]*(AXY + [0.1739; 0.0071; -0.2999]);

%AXY(3,:) = AXY(3,:);% - 9.7932;

%%
% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [q0 q1 q2 q3 B sigma]
x = [0.5, 0.5, 0.5, 0.5]'*ones(1,length(tdgyro));

% Intial P matrix
P = 0.8*eye(length(x(:,1)));
% covariance for w
Q = diag([0.3, 0.3, 0.3, 0.3]);

% covariance for v
RR = [var(AXY(1,:)), var(AXY(2,:)), var(AXY(3,:))]; 

angle = zeros(3,length(tdgyro));

n1 = 1;

for m = 1:1:length(tdgyro)
 
    [y, index, n1] = getyNPVA(tdgyro(m), tdacc, gyro, AXY, n1);
    
    if m == 1
            x(:,m) = x(:,m);
    elseif index == -1
        x(:,m) = x(:,m-1);
    else
        F = getF(gyro,Tgyro, m);

        if m == 1
            x(:,m) = x(:,m);
        else
            x(:,m) = F*x(:,m-1);
        end
        
        P = F*P*F' + Q;
        H = getJacoN(y,x(:,m));

        R = getR(RR);
        % back tracking
        G = P*H'/(H*P*H' + R);%lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m)));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    q0 = x(1,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); q1 = x(2,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); 
    q2 = x(3,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); q3 = x(4,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); 
    
%     angle(1,m) = -atan2(q2*q3-q0*q1, q0*q0+q3*q3-0.5);
%     angle(2,m) = asin(q1*q3+q0*q2);
%     angle(3,m) = -atan2(q1*q2-q0*q3, q0*q0+q1*q1-0.5);
    
    %angle(1,m) = atan2(2*(q2*q3 + q0*q1), q0*q0+q3*q3-q1*q1-q2*q2);
    angle(1,m) = atan2(2*(q2*q3 + q0*q1), q0*q0-q1*q1-q2*q2+q3*q3);%2*(q1*q1+q2*q2)-1);
    %angle(2,m) = asin(-q1*q3+q0*q2);
    angle(2,m) = asin(2*(-q1*q3+q0*q2));
    %angle(3,m) = atan2(q1*q2 + q0*q3, q0*q0+q1*q1-q2*q2-q3*q3);
    angle(3,m) = atan2(2*(q1*q2 + q0*q3), q0*q0+q1*q1-q2*q2-q3*q3);%2*(q2*q2+q3*q3)-1);
    
%     if angle(1,m) > pi/2 || angle(1,m) < -pi/2
%         if angle(2,m) > 0.15
%             angle(2,m) = pi - angle(2,m);
%         elseif angle(2,m) < -0.25
%             angle(2,m) = -pi - angle(2,m);
%         end
%     end

end


angle2 = zeros(3, length(angle)-50);

for m = 1:1:length(tdgyro)-50
    angle2(1,m) = mean(angle(1,m:m+50));
    angle2(2,m) = mean(angle(2,m:m+50));
    angle2(3,m) = mean(angle(3,m:m+50));
    
end

% for m = 2:1:length(angle2)    
% 
% end

figure
subplot(3,1,1), plot(tdgyro(1:length(angle2)), angle2(1,:), 'LineWidth', 2); ylabel('Orientation $\varphi$','Interpreter','latex'); xlabel('time [s]'); title('Estimated orientation along x axis'); grid on; grid minor;
subplot(3,1,2), plot(tdgyro(1:length(angle2)), angle2(2,:), 'LineWidth', 2); ylabel('Orientation $\vartheta$','Interpreter','latex'); xlabel('time [s]');title('Estimated orientation along y axis'); grid on; grid minor;
subplot(3,1,3), plot(tdgyro(1:length(angle2)), angle2(3,:), 'LineWidth', 2); ylabel('Orientation $\psi$','Interpreter','latex'); xlabel('time [s]');title('Estimated orientation along z axis'); grid on; grid minor;
