  function backy=by(x,y,yinput,potinput,pot) 
  D=yinput(1);
  kbt=yinput(2);
  delty=yinput(3);
  deltphi=feval(pot,x,y+delty/2,potinput)-feval(pot,x,y-delty/2,potinput);
  s=deltphi/kbt;
  gam=D/delty^2;
  if abs(s)>1.5e-3
  backy=gam*(-s)/(exp(-s) - 1);
  else
      backy=gam/(1-s*(1/2-s*(1/6-s)));
  end