% Function which solves the shape of an ellipctical flux tube
% in equilibrium for a given initial gradient (=q)
% given by the PDE in "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  Eq. (54) 
%
% input parameters are: x_rho, x_B, A, B1, B2, q, x0

function[zgrit,xequil]=FieldLineEquilibriumStep(xrho,xB,A,B1,B2,q,z0,x0,h,N)

% Integration:
feps=EquilibriumPositionStep(xrho,xB,A,B1,B2,q,x0);

[zgrit,xequil]=RungeKutteSolver(feps,z0, x0, h, N);



    function [dxdz] = EquilibriumPositionStep(xrho,xB,A,B1,B2,q,x0)
        dxdz=@(z,x) A^(-1) *sqrt(A^2*q^2 + ...
            (1+A^2 * q^2)*((x-x0)/((cosh(x0-xrho))^2) ...
            -(tanh(x-xrho)-tanh(x0-xrho)))/( B1-(B2/((cosh(x0-xB))^2)) ) );
    end
end