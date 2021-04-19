clc, clear, close all

sim = 1;

if sim == 1
    
    load('PropGp1DPositionData.mat')
    load('PropGp1DVelocityAccelerationData.mat')

    ta = 0.4657;
    tb = 1.1657;
    
    % Can be used as measurement
    pp = Position_1D(Time >= ta & Time <= tb)'; %-0.99
    vv = MeasVel1D(Time >= ta & Time <= tb)';
    aa = MeasAccel1D(Time >= ta & Time <= tb)';

    td = Time(Time >= ta & Time <= tb)'-ta;
    T = mean(diff(td));

    [PP, VV, AA] = groundtruth1D(td);
    % Start from 0
    PP = PP - PP(1);
    
    % 1D - 3D coordinates
    coordinates = zeros(3, length(PP)*4);

    % x
    coordinates(1,1:length(PP)) = PP;
    coordinates(1,length(PP):length(PP)*3) = PP(end);
    coordinates(1,length(PP)*3+1:length(PP)*4) = flip(PP);

    %y
    coordinates(2,length(PP)+1:length(PP)*2) = PP;
    coordinates(2,length(PP)*2:length(PP)*3) = PP(end);
    coordinates(2,length(PP)*3+1:length(PP)*4) = flip(PP);

    %z
    coordinates(3,length(PP)*2+1:length(PP)*3) = PP;
    coordinates(3,length(PP)*3+1:length(PP)*4) = flip(PP);

    h = axes('Parent', figure());
    hold(h, 'off');


    for i = 1:1:length(PP)*4
        scatter3(h, coordinates(1,i), coordinates(2,i), coordinates(3,i))
        hold on
        plot3(coordinates(1,:), coordinates(2,:), coordinates(3,:))
        xlabel('x(t)')
        ylabel('y(t)')
        zlabel('z(t)')
        grid on
        pause(T)
    end
    
    
else
    
    load('yVectorData.mat')
    
    offset = 0;
    for n = 2:1:length(psif1)
        if psif1(n) - psif1(n-1) < -2
            offset = 2*pi;
        end
        psif1(n) = psif1(n) + offset - pi;  
    end
    
    psif1(1) = psif1(1) - pi;
    
    h = 1e-5;           tc = 0:h:2.00000;
    T = mean(diff(tf)); td = 0:T:2.00000;

    [W_A, W_V, W, XX, rr] = groundtruth2D(tc, td);
    
    % 1D - 3D coordinates
    coordinates = zeros(3, 282);
    
    coordinates(1,1:142) = XX(1,1:142);
    coordinates(1,142+1:282) = XX(1,142);
    
    coordinates(2,:) = XX(4,1:282);
    
    coordinates(3,1:142) = XX(1,142);
    coordinates(3,142+1:282) = XX(1,142+1:282);
    
     
    dt = T; 

    h = axes('Parent', figure());
    hold(h, 'off');


    for i = 1:1:252
        scatter3(h, coordinates(1,i), coordinates(2,i), coordinates(3,i))
        hold on
        plot3(coordinates(1,:), coordinates(2,:), coordinates(3,:))
        pause(dt)
        xlabel('x(t)')
        ylabel('y(t)')
        zlabel('z(t)')
        grid on
    end
    
%     figure
%     subplot(2,1,1),plot(td,XX(1,:),'LineWidth',2), title('Position $x$','interpreter','latex','fontweight','bold'), xlabel('t [s]'),ylabel('x[m]','interpreter','latex'), grid
%     subplot(2,1,2),plot(td,XX(4,:),'LineWidth',2), title('Position $y$','interpreter','latex','fontweight','bold'),xlabel('t [s]'),ylabel('y[m]','interpreter','latex'), grid
%    
    
end