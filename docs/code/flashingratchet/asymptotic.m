function veldeff=asymptotic(param) 
kbt=param(1);          %temperature times boltzmann constant
Dx=param(2);           %diffusion coeff of motor
Dy=param(3);           %diffusion coeff of cargo
yperL=param(4);         %number of y mesh points in a period
xmax=param(5);          %number of x values in a period
periodtotal=param(6);   %number of periods y goes through
periodsahead=param(7);  %number of periods the mesh for y is calculated beyond origin (in positive direction) 
K=param(8);             %spring constant
F=param(9);              %positive forcing on motor
A=param(10);             %values for amplitude 
period=param(11);        %length of period
k1=param(12);            %switching rate from potential one to potential two
k2=param(13);	         %switching rate from potential two to potential one
pot1=@phia ;             %function "phia()" will be potential one--allows us to pass the function to another function
pot2=@phi ;              %function "phi()" will be potential two

%setup vector to pass to potential
potinput=zeros(1,4);
potinput(1)=K;
potinput(2)=F;
potinput(3)=A;
potinput(4)=period;

%create y mesh
N=yperL*periodtotal; %total number of mesh points for y
y=zeros(1,N);
delty=period/yperL;  %delta y for the ymesh 
for i=(N+2-periodsahead*yperL):N
    y(i)=y(i-1)+delty;
end %for i
for i=(N-periodsahead*yperL):-1:1
    y(i)=y(i+1)-delty;
end %for i

%create yinput vector
yinput=zeros(1,3);
yinput(1)=Dy;
yinput(2)=kbt;
yinput(3)=delty;


%create x mesh
x=zeros(1,xmax);
deltx=period/xmax; %delta x for the x mesh
for i=1:xmax-1
    x(i+1)=x(i)+deltx;
end %for i

%create xinput vector
xinput=zeros(1,3);
xinput(1)=Dx;
xinput(2)=kbt;
xinput(3)=deltx;


%Lprime
Lprime1=zeros(xmax*N);
for n=1:N
for i=1:xmax
    if i<xmax
    Lprime1(i+(n-1)*xmax,i+1+(n-1)*xmax)=fx(x(i),y(n),xinput,potinput,pot1);
    end %if
    if i>1
    Lprime1(i+(n-1)*xmax,i+(n-1)*xmax-1)=bx(x(i),y(n),xinput,potinput,pot1) ;
    end %if
end % for i
for i=1:xmax    
   if n<N
   Lprime1(i+n*xmax,i+(n-1)*xmax)=by(x(i),y(n+1),yinput,potinput,pot1); %submatrix (n+1,n)
   Lprime1(i+(n-1)*xmax,i+n*xmax)=fy(x(i),y(n),yinput,potinput,pot1);   %submatrix (n,n+1)
   end %if
end %for i
end % for n


%Set up L+'
Lplusprime1=zeros(xmax*N);
for n=1:yperL
    Lplusprime1(n*xmax,1)=fx(x(xmax),y(n),xinput,potinput,pot1);
end %for n
for n=(yperL+1):N
    Lplusprime1(n*xmax,xmax*(n-yperL)-(xmax-1))=fx(x(xmax),y(n),xinput,potinput,pot1);
end %for n

%Set up L-'
Lminusprime1=zeros(xmax*N);
for n=N:-1:(N-yperL+1)
    Lminusprime1(n*xmax-(xmax-1),xmax*N)=bx(x(1),y(n),xinput,potinput,pot1);
end %for n
for n=(N-yperL):-1:1
    Lminusprime1(n*xmax-(xmax-1),xmax*(n+yperL))=bx(x(1),y(n),xinput,potinput,pot1);
end %for n


%Lprime
Lprime2=zeros(xmax*N);
for n=1:N
for i=1:xmax
    if i<xmax
    Lprime2(i+(n-1)*xmax,i+1+(n-1)*xmax)=fx(x(i),y(n),xinput,potinput,pot2);
    end %if
    if i>1
    Lprime2(i+(n-1)*xmax,i+(n-1)*xmax-1)=bx(x(i),y(n),xinput,potinput,pot2) ;
    end %if
end % for i
for i=1:xmax    
   if n<N
   Lprime2(i+n*xmax,i+(n-1)*xmax)=by(x(i),y(n+1),yinput,potinput,pot2); %submatrix (n+1,n)
   Lprime2(i+(n-1)*xmax,i+n*xmax)=fy(x(i),y(n),yinput,potinput,pot2);   %submatrix (n,n+1)
   end %if
end %for i
end % for n n


%Set up L+'
Lplusprime2=zeros(xmax*N);
for n=1:yperL
    Lplusprime2(n*xmax,1)=fx(x(xmax),y(n),xinput,potinput,pot2);
end %for n
for n=(yperL+1):N
    Lplusprime2(n*xmax,xmax*(n-yperL)-(xmax-1))=fx(x(xmax),y(n),xinput,potinput,pot2);
end %for n

%Set up L-'
Lminusprime2=zeros(xmax*N);
for n=N:-1:(N-yperL+1)
    Lminusprime2(n*xmax-(xmax-1),xmax*N)=bx(x(xmax-(xmax-1)),y(n),xinput,potinput,pot2);
end %for n
for n=(N-yperL):-1:1
    Lminusprime2(n*xmax-(xmax-1),xmax*(n+yperL))=bx(x(xmax-(xmax-1)),y(n),xinput,potinput,pot2);
end %for n


for i=1:N
    for j=1:xmax
       %s=feval(potx1,x(j),y(i),xinput)-feval(poty1,x(j),y(i),K);
        k12((i-1)*xmax+j,(i-1)*xmax+j)=k1; %*exp(s);
    end % for j
end %for i


for i=1:N
    for j=1:xmax
        %s=feval(potx2,x(j),y(i),xinput)-feval(poty2,x(j),y(i),K);
        k21((i-1)*xmax+j,(i-1)*xmax+j)=k2; %*exp(s);
    end %for j
end %for i


Lprime=[Lprime1,k12;k21,Lprime2];
Lplusprime=[Lplusprime1,zeros(length(Lplusprime1));zeros(length(Lplusprime1)),Lplusprime2];
Lminusprime=[Lminusprime1,zeros(length(Lminusprime1));zeros(length(Lminusprime1)),Lminusprime2];


%make diagonal for L matrix
for i=1:length(Lprime)
    Lprime(i,i)=-sum(Lprime(i,:))-sum(Lplusprime(i,:))-sum(Lminusprime(i,:));
end %for i

L=Lprime'; %transpose all of matrices
Lplus=Lplusprime';
Lminus=Lminusprime';
M=L+Lplus+Lminus;
Mtemp=M;
Mtemp(2*N*xmax,:)=ones(1,2*N*xmax); %replace last row of M with ones (add condition)
rhs=zeros(2*N*xmax,1);
rhs(2*N*xmax)=1; %condition for ps
s=sparse(Mtemp);
ps=s\rhs;

%asymptotic velocity
velocity=period*sum((Lplus-Lminus)*ps);
rhs=sum((Lplus-Lminus)*ps)*ps-(Lplus-Lminus)*ps;
rhs(2*N*xmax)=0; %condition for r
r=s\rhs;

%effective diffusion
deff=0.5*period^2*sum((Lplus+Lminus)*ps+2*(Lplus-Lminus)*r);
veldeff=[velocity,deff];
