% Algorithm which solves the shape of an ellipctical flux tube
% in equilibrium 
% given by the PDE in "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  Eq. (54) 
%
% input parameters are: x_rho, x_B, A, B1, B2, q, x0
clear all
close all

% input parameters:
A = 0.161604;
%x0 = 1.2;
%x0=1.1;
x0=1.0;
%x0=0.5;
xrho = 2;
xB = 0.8;
B1 = 3;
B2 = 1.8;
%q=3.25;
%q=4.5;
q=6.15;
%q = 25.25;

h=0.05;
zmax=0.5;
N=zmax/h;
z0=0;

[zgritQ,xequilQ]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,h,N);
plot(zgritQ,xequilQ)




