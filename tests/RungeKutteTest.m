% Code which tests RungeKutteSolver (Or other 1st order ODE solvers)
% Calculates solution of RungeKutteSolver with some analytic solution
addpath('../tools');
clear all
close all
f=@(x,y) x+y;
x0=0;
y0=0;
h=0.2;
xmax=1;
N=xmax/h;

[xgrit,yRK]=RungeKutteSolver(f,x0,y0,h,N);

yCorrectValues=exp(xgrit)-xgrit-1;

plot(xgrit,yRK,'--',xgrit,yCorrectValues)

