function backx=bx(x,y,xinput,potinput) 
  D=xinput(1);
  kbt=xinput(2);
  deltx=xinput(3);
  deltphi=phi(x+deltx/2,y,potinput)-phi(x-deltx/2,y,potinput);
  s=deltphi/kbt;
  gam=D/deltx^2;
  if abs(s)>1.5e-3
  backx=gam*(-s)/(exp(-s) - 1);
  else
      backx=gam/(1-s*(1/2-s*(1/6-s)));
  end
 
  
  
 