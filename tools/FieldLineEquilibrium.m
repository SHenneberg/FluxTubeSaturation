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
q=4.23;
%q = 25.25;

h=0.01;
zmax=1.; %DON'T CHANGE this anymore!! It's important for the root finding algorithm
N=zmax/h;
z0=0; %DON'T CHANGE this anymore!! It's important for the root finding algorithm
%input for root finding
qinit1=q+1;
qinit2=q-1;
itr=20;
tol=10^(-8);


xsym=@(qvar) xsymmetry(xrho,xB,A,B1,B2,z0,qvar,x0,h,N);


[qroot,success,ifinal] = SecantMethod(xsym,qinit1,qinit2,itr,tol);

qroot
ifinal

% second order ODE:
[zgritQ,xequilQ]=FieldLineEquilibriumStep2(xrho,xB,A,B1,B2,qroot,z0,x0,h,N);
plot(zgritQ,xequilQ)

% first order ODE:
%[zgritQ1,xequilQ1]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,h,N);
%plot(zgritQ1,xequilQ1,'--',zgritQ,xequilQ)

% function which need to be 0 when the right q is found:
% symmetry of x with respect to z=0.5
function[xsym]=xsymmetry(xrho,xB,A,B1,B2,z0,q,x0,h,N)
    [zgritQ,xsymvar]=FieldLineEquilibriumStep2(xrho,xB,A,B1,B2,q,z0,x0,h,N);
    xsym=xsymvar(1)-xsymvar(N+1);    
end


