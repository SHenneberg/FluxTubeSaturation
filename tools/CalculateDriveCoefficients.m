% Calculate Drive coefficients

clear all
close all

% input parameters:
A = 0.161604;
%AB1=0.07834; % Value in the paper
%AB2=0.04701; % Value in the paper
AB1=3.;
AB2=2.972;
x0min=0.2;
x0max=2.2;
x0step=0.05;
%xrho = 2;
xrho =1.1;
xB = 0.8;
%B1 = 3; % Value in the paper
B1=AB1/A^2
%B2 = 1.8; % Value in the paper
B2=AB2/A^2
q=4.23;

x0scan=x0min:x0step:x0max;

x0=0.8
[gamma,C4]=DriveCoefficients(A,B1, B2, x0, xB, xrho)
[gammaS,C4S]=DriveCoefficients(A,B1, B2, x0scan, xB, xrho);

plot(x0scan,gammaS, '--',x0scan,C4S)