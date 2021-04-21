clc, clear, close all

%load coord3.mat;

% cm=1;
% in=0.393*cm
% pr=in*60000%*-1

s=serial('/dev/tty.usbserial-1440','BaudRate',19200,'Parity', 'none', 'DataBits',8, 'StopBits', 1, 'terminator', 'CR');
fopen(s);
%pause(2);

%fprintf(s,'PRX=1000')

% x: -32000 - 190,000
% y: -14,0000 - 1000
% z:    0 - 24,000

x = ones(7,3);        
x(1,1) = 0;       x(1,2) = 0;      x(1,3) = 0;
x(2,1) = 10000;   x(2,2) = 0;      x(2,3) = 0;
x(3,1) = 0;       x(3,2) = 0;      x(3,3) = 0;

x(4,1) = 120000;  x(4,2) = 0;      x(4,3) = 0;
x(5,1) = 120000;  x(5,2) = 92000; x(5,3) = 0;
x(6,1) = 120000;  x(6,2) = 92000; x(6,3) = 14000;
x(7,1) = 0;       x(7,2) = 0;      x(7,3) = 0;

% xdot = ones(3,1);  xdot(1) = 10000; xdot(2) = 40000; xdot(3) = 4000;
% 
% y = ones(3,1);        y(1) = -5000; y(2) = 140000; y(3) = -3000;
% 
% z = ones(3,1);        z(1) = -5000; z(2) = 140000; z(3) = 20000;
% x: -180000 - -3,000
% y:   0     -140,000
% z:  -3,000 - 17,000

% -180000,0,-3000

%% start
fprintf(s,'SH X')
i = 2; %Initialize counter
n=6; 
%i = i
%dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
%fprintf(s,ac);
%fprintf(s,dc);
fprintf(s,'BG X');%BEGIN Motion along y axis only
pause(5);
i=i+1;

%% back
fprintf(s,'SH X')
%i = i
%dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
%fprintf(s,ac);
%fprintf(s,dc);
fprintf(s,'BG X');%BEGIN Motion along y axis only
pause(5);
i=i+1;


%% X
fprintf(s,'SH X')
%i = i
%dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG X');%BEGIN Motion along y axis only
pause(25);
i=i+1;

%% Y
s=serial('/dev/cu.usbserial-1460','BaudRate',19200,'Parity', 'none', 'DataBits',8, 'StopBits', 1, 'terminator', 'CR');
fopen(s);

fprintf(s,'SH Y')

i = i
%dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG Y');%BEGIN Motion along y axis only
pause(25);
i=i+1;

%% Z
fprintf(s,'SH Z')
 
%i = i
%dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG Z');%BEGIN Motion along y axis only
pause(15);
i=i+1;

%% XYZ
% fprintf(s,'SH XYZ')
%  
% %i = i
% %dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
% dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
% 
% sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
% ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
% dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
% 
% %re = sprintf('TP Y') ;
% 
% fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
% fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
% fprintf(s,ac);
% fprintf(s,dc);
% fprintf(s,'BG XYZ');%BEGIN Motion along y axis only
% pause(20);


%i = 0; %Initialize counter
% pause(0.5);
% while i<n, 
%     dis=sprintf('PR %i,%i,%i',pr,0,0)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%     sp=sprintf('SP %i,%i,%i',5000,0,0)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%     fprintf(s,dis);%WRITE DISTANCE string to Controller through serial port 
%     fprintf(s,sp);%WRITE DISTANCE string to Controller through serial port 
%     fprintf(s,'BG X');%BEGIN Motion along y axis only
%     pause(0.5); 
%     i=i+1;
% end 
fclose(s)