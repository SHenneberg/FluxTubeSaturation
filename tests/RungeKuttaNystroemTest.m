% Code which tests RungeKutteNystroemSolver (Or other 2nd order ODE solvers)
% Calculates solution of RungeKutteNystroemSolver with some analytic solution
addpath('../tools');
clear all
close all
f=@(x,y,dy) x+y;
x0=0;
y0=1;
dy0=-2;
h=0.05;
xmax=15.;
N=xmax/h;

[xgrit,yRK]=RungeKuttaNystroemSolver(f,x0,y0,dy0,h,N);

for yloop=1:N+1
    yCorrectValues(yloop)=-exp(-xgrit(yloop))*(-1+xgrit(yloop)*exp(xgrit(yloop)));
end


plot(xgrit,yRK,'--',xgrit,yCorrectValues)