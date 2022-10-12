clear all;
format long;
param=zeros(1,11);
param(1)=4.2;         %temperature times boltzmann constant
param(2)=9e3 ;        %diffusion coeff of motor
param(3)=9e2;         %diffusion coeff of cargo
param(4)=1;           %number of y mesh points in a period
param(5)=50;          %number of x values in a period
param(6)=3;           %number of periods y goes through
param(7)=2;  	      %number of periods the mesh for y is calculated ahead of the origin (in positive direction) 
param(8)=0; 	      %spring constant K 
param(9)=10;          %positive forcing on motor
param(10)=8*param(1);           %values for amplitude 
param(11)=8;     		%period for the motor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is how "asymptotic()" uses the parameter values.
% kbt=param(1);          %temperature times boltzmann constant
% Dx=param(2);           %diffusion coeff of motor
% Dy=param(3);           %diffusion coeff of cargo
% yperL=param(4);         %number of y mesh points in a period
% xmax=param(5);          %number of x values in a period
% periodtotal=param(6);   %number of periods y goes through
% periodsahead=param(7);  %number of periods the mesh for y is calculated beyond origin (in positive direction) 
% K=param(8);             %spring constant
% F=param(9);              %positive forcing on motor
% A=param(10);             %values for amplitude 
% period=param(11);        %length of period

returnvec=asymptotic(param);
   velocity=returnvec(1);
   deff=returnvec(2);
