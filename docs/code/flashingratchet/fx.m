 function forwardx=fx(x,y,xinput,potinput,pot) 
  D=xinput(1);
  kbt=xinput(2);
  deltx=xinput(3);
  deltphi=feval(pot,x+deltx/2,y,potinput)-feval(pot,x-deltx/2,y,potinput);
  s=deltphi/kbt;
  gam=D/deltx^2;
  if abs(s) > 1.5e-3
      forwardx=gam*s/(exp(s)-1);
  else
      forwardx=gam/(s*(s*(s/24+1/6)+1/2)+1);
  end
  