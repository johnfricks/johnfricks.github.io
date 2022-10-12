function pot=phi(x,y,potinput)
    K=potinput(1);
    F=potinput(2);
   pot=K/2*(y-x)^2+F*y;
