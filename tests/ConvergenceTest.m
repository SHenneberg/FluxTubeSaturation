% Convergence test
% input: f=solver (ode solver) + dx=stepsize + h=initial step + Nconv=number of steps
addpath('../tools');

%input parameters for ode solver:
A = 0.161604;
x0 = 0.5;
xrho = 2;
xB = 0.8;
B1 = 3;
B2 = 1.8;
q = 5;
h=0.01;
zmax=0.5;
N=zmax/h;
z0=0;

%f: [zgrit,xequil]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,dx,N);

%step size dx;
dx=0.01;
%initial finite difference: h=0.05
% number of steps:
Nconv=5;


for k=1:Nconv
    hconv(k)=h+dx*(k-1);
    [zgrit,xequil]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,hconv(k),N);
    xequilCONV(k)=xequil(N+1);
end

plot(hconv,xequilCONV)


