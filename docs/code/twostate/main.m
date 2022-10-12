clear all;
format long;
kbt=4.2;         %temperature times boltzman
Dy=0.001*3e5;    %diffusion coeff of cargo            
period=8;        %length of period
yperL=10;        %number of y mesh points in a period
xmax=2;          %number of x values in L
periodtotal=4;   %number of periods y goes through
periodsahead=2;  %number of periods the mesh for y is calculated beyond origin (in positive direction) 
K=1;             %array of values for K 
amplitude=0;     %array of values for A 
atp=20; 	 %ATP concentration

f1=3.75*atp;	%forward step rate of motor alone to move from state 1 to state 2
f2=141.1;	%forward step rate of motor alone to move from state 2 to state 1 (of the next period)
b1=0.034;	%backward step rate of motor alone to move from state 1 to state 2 (of previous period)
b2=0.0047;      %backward step rate of motor alone to move from state 2 to state 1 
F=0;   		%backward force felt by the bead(the cargo)

%the following sets up a parameter vector to pass to the function "asymptotic()" 
param=zeros(1,13); 
param(1)=kbt;
param(2)=b1;
param(3)=b2;
param(4)=f1;
param(5)=f2;
param(6)=Dy;
param(7)=yperL;
param(8)=xmax;
param(9)=periodtotal;
param(10)=periodsahead;
param(11)=K;
param(12)=F;
param(13)=period;


%call function "asymptotic()" which returns a vector of length 2 with the asymptotic velocity and the effective diffusion
returnvec=asymptotic(param);
vel=returnvec(1);
deff=returnvec(2);



