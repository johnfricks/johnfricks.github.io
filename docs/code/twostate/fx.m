function forwardx=fx(x,y,xinput,potinput) 
  kbt=xinput(1);
  deltx=xinput(2);
  f1=xinput(5);
  f2=xinput(6);
  if  x< 1e-3 
      forwardx=f1*exp(0.5/kbt*(phi(x,y,potinput)-phi(x+deltx,y,potinput)));
  else
      forwardx=f2*exp(0.13/kbt*(phi(x,y,potinput)-phi(x+deltx,y,potinput)));
  end
