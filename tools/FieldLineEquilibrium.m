% Algorithm which solves the shape of an ellipctical flux tube
% in equilibrium 
% given by the PDE in "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  Eq. (54) 
%
% input parameters are: x_rho, x_B, A, B1, B2, q, x0

% default input parameters:
A = 0.161604;
x0 = 0.5;
xrho = 2;
xB = 0.8;
B1 = 3;
B2 = 1.8;
q = 3.25;

h=0.005;
zmax=0.5;
N=zmax/h;
z0=0;

[zgrit,xequil]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,h,N);
plot(zgrit,xequil)




