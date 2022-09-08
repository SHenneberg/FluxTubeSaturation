% Algorithm which solves the shape of an ellipctical flux tube
% in equilibrium 
% given by the PDE in "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  Eq. (54) 
%
% input parameters are: x_rho, x_B, A, B1, B2, q, x0

% default input parameters:
A = 0.161604;
x0 = 1.2;
xrho = 2;
xB = 0.8;
B1 = 3;
B2 = 1.8;
q = 3.25;

h=0.01;
zmax=1;
N=zmax/h;
z0=0

% Integration:
f=EquilibriumPositionStep(xrho,xB,A,B1,B2,q,x0);
[zgrit,xequil]=RungeKutteSolver(f,z0, x0, h, N);

plot(zgrit,xequil)




function [dxdz] = EquilibriumPositionStep(xrho,xB,A,B1,B2,q,x0)
    dxdz=@(z,x) A^(-1) *sqrt(A^2*q^2 + (1+A^2 * q^2)*((x-x0)/(cosh(x0-xrho))  - (tanh(x-xrho)-tanh(x0-xrho))/(B1-(B2/(cosh(x0-xB)^2)))))
end

