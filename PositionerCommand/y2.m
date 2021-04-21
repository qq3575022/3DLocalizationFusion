s=serial('/dev/tty.usbserial-1440','BaudRate',19200,'Parity', 'none', 'DataBits',8, 'StopBits', 1, 'terminator', 'CR');
fopen(s);

x = ones(7,3);        
x(1,1) = 0;       x(1,2) = 0;      x(1,3) = 0;
x(2,1) = 0;       x(2,2) = 2000;      x(2,3) = 0;
x(3,1) = 0;       x(3,2) = 0;      x(3,3) = 0;

x(4,1) = 120000;  x(4,2) = 0;      x(4,3) = 0;
x(5,1) = 120000;  x(5,2) = 92000; x(5,3) = 0;
x(6,1) = 120000;  x(6,2) = 92000; x(6,3) = 14000;
x(7,1) = 0;       x(7,2) = 0;      x(7,3) = 0;

i = 2;

fprintf(s,'SH Y')

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
fprintf(s,'BG Y');%BEGIN Motion along y axis only
pause(5);
