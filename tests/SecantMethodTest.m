% Root Finding test for secant method

addpath('../tools');
clear all
close all

%function
f=@(x) x.^2-2;
%input for root finding
x0=1;
x1=1.2;
itr=10;
tol=10^(-10);

[root,success,itr] = SecantMethod(f,x0,x1,itr,tol)

xgrit=linspace(-3,3);
plot(xgrit,f(xgrit))


