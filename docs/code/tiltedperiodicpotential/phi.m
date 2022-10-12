function pot=phi(x,y,input) 
K=input(1);
F=input(2);
A=input(3);
L=input(4);
pot=A*sin(2*pi/L*x)-F*x+K/2*(y-x)^2;