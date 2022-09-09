% Function which solves the shape of an ellipctical flux tube
% in equilibrium for a given initial gradient (=q)
% given by the 2nd order PDE in "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  derived from Eq. (34) with normalisations
% from (48)-(52)
%
% input parameters are: x_rho, x_B, A, B1, B2, q, x0

% first order pde 
%2nd order PDE
function[zgrit,xequil]=FieldLineEquilibriumStep2(xrho,xB,A,B1,B2,q,z0,x0,h,N)

% Integration:
feps=EquilibriumPositionStep2(xrho,xB,A,B1,B2,x0);

[zgrit,xequil,dxequil]=RungeKuttaNystroemSolver(feps,z0, x0,q, h, N);



    function [d2xdz2] = EquilibriumPositionStep2(xrho,xB,A,B1,B2,x0)
        d2xdz2=@(z,x,dx) 0.5*A^(-2) *...
            (1+A^2*dx^2)/(B1-B2/((cosh(x0-xB))^2)+(x-x0)/((cosh(x0-xrho))^2) ...
            - (tanh(x-xrho)-tanh(x0-xrho)))^2  * (sech(x0-xrho)^2-sech(x-xrho)^2);
    end
end
