clc, clear, close all

[mag,gyro,AXY,Tmag,tdmag,Tgyro,tdgyro,Tacc,tdacc,T3,td3] = meas3D();


% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]
% State vector [q0 q1 q2 q3 B sigma]
qx = [0.5, 0.5, 0.5, 0.5]'*ones(1,length(tdgyro));

% Intial P matrix
qP = 0.8*eye(length(qx(:,1)));
% covariance for w
qQ = diag([0.3, 0.3, 0.3, 0.3]);

% covariance for v
qRR = [var(AXY(1,:)), var(AXY(2,:)), var(AXY(3,:))]; 

angle = zeros(3,length(tdgyro));

n1 = 1;

% for m = 1:1:length(tdgyro)
%  
%     [qy, index, n1] = getqyNPVA(tdgyro(m), tdacc, gyro, AXY, n1);
% 
%     if m == 1
%             qx(:,m) = qx(:,m);
%     elseif index == -1
%         qx(:,m) = qx(:,m-1);
%     else
%         qF = getqF(gyro,Tgyro, m);
% 
%         qx(:,m) = qF*qx(:,m-1);
% 
%         qP = qF*qP*qF' + qQ;
%         qH = getJacoqN(qy,qx(:,m));
% 
%         qR = getqR(qRR);
%         % back tracking
%         qG = qP*qH'/(qH*qP*qH' + qR);%lambda
%         % e
%         qx(:,m) = qx(:,m) + qG*(qy - getqH(qx(:,m)));
%         qP = (eye(length(qx(:,1))) - qG*qH)*qP;
%     end
%     
%     q0 = qx(1,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); q1 = qx(2,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); 
%     q2 = qx(3,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); q3 = qx(4,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); 
%     
% %     angle(1,m) = -atan2(q2*q3-q0*q1, q0*q0+q3*q3-0.5);
% %     angle(2,m) = asin(q1*q3+q0*q2);
% %     angle(3,m) = -atan2(q1*q2-q0*q3, q0*q0+q1*q1-0.5);
%     
%     angle(1,m) = atan2(-2*(q2*q3 + q0*q1), q0*q0-q1*q1-q2*q2+q3*q3);%2*(q1*q1+q2*q2)-1);
%     angle(2,m) = asin(-2*(q1*q3-q0*q2));
%     angle(3,m) = atan2(-2*(q1*q2 + q0*q3), q0*q0+q1*q1-q2*q2-q3*q3);%2*(q2*q2+q3*q3)-1);
%     
% end
% 
% angle2 = zeros(3, length(angle)-50);
% 
% for m = 1:1:length(tdgyro)-50
%     angle2(1,m) = mean(angle(1,m:m+50));
%     angle2(2,m) = mean(angle(2,m:m+50));
%     angle2(3,m) = mean(angle(3,m:m+50));
%     
% end

% for m = 2:1:length(angle2)    
% 
% end
% figure
% subplot(311), plot(x(1,:))
% subplot(312), plot(x(2,:))
% subplot(313), plot(x(3,:))

% figure
% subplot(3,1,1), plot(tdgyro(1:length(angle2)), angle2(1,:), 'LineWidth', 2); ylabel('Orientation $\varphi$','Interpreter','latex'); xlabel('time [s]'); title('Estimated orientation along x axis'); grid on; grid minor;
% subplot(3,1,2), plot(tdgyro(1:length(angle2)), angle2(2,:), 'LineWidth', 2); ylabel('Orientation $\vartheta$','Interpreter','latex'); xlabel('time [s]');title('Estimated orientation along y axis'); grid on; grid minor;
% subplot(3,1,3), plot(tdgyro(1:length(angle2)), angle2(3,:), 'LineWidth', 2); ylabel('Orientation $\psi$','Interpreter','latex'); xlabel('time [s]');title('Estimated orientation along z axis'); grid on; grid minor;

%%

% State vector [xdotdot, ydotdot, orientation, angular velocity]
x = [0.9944; 0.0623; 0.0054; 0.9999; -0.004; 0.9880; 0.1717; 0.008; -0.3]*ones(1,length(tdgyro));

% Intial P matrix
P = 0.8*eye(length(x(:,1)));
% covariance for w
Q = diag([0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]);


% covariance for v
%RR = [0.01^2, 0.01^2, 0.01^2, 0.2^2, 0.014^2, 0.2^2, 0.014^2, 0.2^2, 0.014^2]; 
RR = [0.8, 0.8, 0.8]; % can tune

for m = 1:1:length(tdgyro)
    
    [qy, index, n1] = getqyNPVA(tdgyro(m), tdacc, gyro, AXY, n1);
     y = getyNPVA(m, AXY);
     
    if m == 1
        
            qx(:,m) = qx(:,m);
            x(:,m) = x(:,m);
            
    elseif index == -1
        qx(:,m) = qx(:,m-1);
    else
        qF = getqF(gyro,Tgyro, m);

        qx(:,m) = qF*qx(:,m-1);

        qP = qF*qP*qF' + qQ;
        qH = getJacoqN(qy,qx(:,m));

        qR = getqR(qRR);
        % back tracking
        qG = qP*qH'/(qH*qP*qH' + qR);%lambda
        % e
        qx(:,m) = qx(:,m) + qG*(qy - getqH(qx(:,m)));
        qP = (eye(length(qx(:,1))) - qG*qH)*qP;
        
        %
        q0 = qx(1,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); q1 = qx(2,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); 
        q2 = qx(3,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); q3 = qx(4,m)/sqrt(qx(1,m)^2+qx(2,m)^2+qx(3,m)^2+qx(4,m)^2); 
        %
        
        x(:,m) = x(:,m-1);

        P = P + Q;
        a = [2*(q1*q3 - q0*q2); 2*(q2*q3 + q0*q1); q0*q0 - q1*q1 - q2*q2 + q3*q3]*9.793;
        
        H = getJacoN(y,x(:,m),a);

        % back tracking
        G = P*H'/(H*P*H' + RR);%lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m), a));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    
end
biasA = [mean(x(1,:)), mean(x(2,:)), mean(x(3,:)), mean(x(4,:)), mean(x(5,:)), mean(x(6,:)), mean(x(7,:)), mean(x(8,:)), mean(x(9,:))]
% biasG = [mean(x(11,:)), mean(x(13,:)), mean(x(15,:))]
% biasM = [mean(x(10,:)), mean(x(12,:)), mean(x(14,:))]
% 
% driftM = [mean(x(16,:)), mean(x(17,:)), mean(x(18,:))]
%%
figure
subplot(9,1,1), plot(tdgyro, x(1,:), 'b', 'LineWidth', 2); ylim([-0.1, 1]);%hold on; plot(t_x, ones(1, length(x(1,:)))*mean(x(1,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean'); title('Estimated bias of acceleration along x axis'); grid on;
subplot(9,1,2), plot(tdgyro, x(2,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.1]);%hold on; plot(t_x, ones(1, length(x(2,:)))*mean(x(2,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean'); title('Estimated bias of acceleration along y axis'); grid on;
subplot(9,1,3), plot(tdgyro, x(3,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.1]);%hold on; plot(t_x, ones(1, length(x(3,:)))*mean(x(3,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean'); title('Estimated bias of acceleration along z axis'); grid on;
subplot(9,1,4), plot(tdgyro, x(4,:), 'b', 'LineWidth', 2); ylim([-0.1, 1]);%hold on;plot(t_x, ones(1, length(x(4,:)))*mean(x(4,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of magnetometor along x axis'); grid on;
subplot(9,1,5), plot(tdgyro, x(5,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.2]);%hold on;plot(t_x, ones(1, length(x(5,:)))*mean(x(5,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of magnetometor along y axis'); grid on;
subplot(9,1,6), plot(tdgyro, x(6,:), 'b', 'LineWidth', 2); ylim([-0.1, 1]);%hold on;plot(t_x, ones(1, length(x(6,:)))*mean(x(6,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of magnetometor along z axis'); grid on;
subplot(9,1,7), plot(tdgyro, x(7,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.1]);%hold on;plot(t_x, ones(1, length(x(7,:)))*mean(x(7,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of gyroscope along x axis'); grid on;
subplot(9,1,8), plot(tdgyro, x(8,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.1]);%hold on;plot(t_x, ones(1, length(x(8,:)))*mean(x(8,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of gyroscope along y axis'); grid on;
subplot(9,1,9), plot(tdgyro, x(9,:), 'b', 'LineWidth', 2); ylim([-0.1, 0.3]);%hold on;plot(t_x, ones(1, length(x(9,:)))*mean(x(9,:)), 'r', 'LineWidth', 2); legend('estimated', 'mean');title('Estimated bias of gyroscope along z axis'); grid on;