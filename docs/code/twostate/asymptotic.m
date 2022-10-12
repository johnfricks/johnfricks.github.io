function veldeff=asymptotic(param) 
kbt=param(1);          %temperature times boltzmann constant
b1=param(2);           %backward step rate of motor alone to move from state 1 to state 2 (of previous period)
b2=param(3);           %backward step rate of motor alone to move from state 2 to state 1 
f1=param(4);           %forward step rate of motor alone to move from state 1 to state 2
f2=param(5);		%forward step rate of motor alone to move from state 2 to state 1 (of the next period)
Dy=param(6);           %diffusion coeff of cargo
yperL=param(7);         %number of y mesh points in a period
xmax=param(8);          %number of x values in a period
periodtotal=param(9);   %number of periods y goes through
periodsahead=param(10);  %number of periods the mesh for y is calculated beyond origin (in positive direction) 
K=param(11);             %spring constant
F=param(12);              %negative forcing on bead
period=param(13);        %length of period


%setup vector to pass to potential
potinput=zeros(1,2);
potinput(1)=K;
potinput(2)=F;

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
xinput=zeros(1,6);
xinput(1)=kbt;
xinput(2)=deltx;
xinput(3)=b1;
xinput(4)=b2;
xinput(5)=f1;
xinput(6)=f2;


%set up L'
Lprime=zeros(xmax*N);
for n=1:N
for i=1:xmax
    if i<xmax
    Lprime(i+(n-1)*xmax,i+1+(n-1)*xmax)=fx(x(i),y(n),xinput,potinput);
    end %if
    if i>1
    Lprime(i+(n-1)*xmax,i+(n-1)*xmax-1)=bx(x(i),y(n),xinput,potinput);
    end %if
end % for i
for i=1:xmax    
   if n<N
   Lprime(i+n*xmax,i+(n-1)*xmax)=by(x(i),y(n+1),yinput,potinput); %submatrix (n+1,n)
   Lprime(i+(n-1)*xmax,i+n*xmax)=fy(x(i),y(n),yinput,potinput);   %submatrix (n,n+1)
   end
end %for i
end %for n

%Set up L+'
Lplusprime=zeros(xmax*N);
for n=1:yperL
    Lplusprime(n*xmax,1)=fx(x(xmax),y(n),xinput,potinput);
end %for n
for n=(yperL+1):N
    Lplusprime(n*xmax,xmax*(n-yperL)-(xmax-1))=fx(x(xmax),y(n),xinput,potinput);
end %for n

%Set up L-'
Lminusprime=zeros(xmax*N);
for n=N:-1:(N-yperL+1)
    Lminusprime(n*xmax-(xmax-1),xmax*N)=bx(x(1),y(n),xinput,potinput);
end %for n
for n=(N-yperL):-1:1
    Lminusprime(n*xmax-(xmax-1),xmax*(n+yperL))=bx(x(1),y(n),xinput,potinput);
end %for n

%Calculate diagonal of L'
for i=1:(N*xmax)
    Lprime(i,i)=-sum(Lprime(i,:))-sum(Lplusprime(i,:))-sum(Lminusprime(i,:));
end %for i

%transpose all of matrices
L=Lprime'; 
Lplus=Lplusprime';
Lminus=Lminusprime';

M=L+Lplus+Lminus;
Mtemp=M;
Mtemp(N*xmax,:)=ones(1,N*xmax); %replace last row of M with ones (add condition)
rhs=zeros(N*xmax,1);
rhs(N*xmax)=1; %condition for ps
ps=Mtemp\rhs;
%asymptotic velocity
velocity=period*sum((Lplus-Lminus)*ps);

rhs=sum((Lplus-Lminus)*ps)*ps-(Lplus-Lminus)*ps;
rhs(N*xmax)=0; %condition for r
r=Mtemp\rhs;

%effective diffusion
deff=0.5*period^2*sum((Lplus+Lminus)*ps+2*(Lplus-Lminus)*r);
veldeff=[velocity, deff];