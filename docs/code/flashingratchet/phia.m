function pot=phia(x,y,potinput) 
K=potinput(1);
F=potinput(2);
A=potinput(3);
L=potinput(4);
xtemp=-x;
pot=K/2*(y-x)^2+A/(2*pi)*(sin(2*pi/L*xtemp)-1/2*sin(4*pi/L*xtemp)+1/3*sin(6*pi/L*xtemp))-F*xtemp;
