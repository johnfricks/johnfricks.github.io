function backx=bx(x,y,xinput,potinput) 
  kbt=xinput(1);
  deltx=xinput(2);
  b1=xinput(3);
  b2=xinput(4);
  if x<1.5e-3
  backx= b1*exp(0.87/kbt*(phi(x,y,potinput)-phi(x-deltx,y,potinput))) ;
  else
  backx= b2*exp(0.5/kbt*(phi(x,y,potinput)-phi(x-deltx,y,potinput))) ;
  end
 
