% Function + code which calculates the linear and nonlinear drive coefficients
% for a specific normalised equilibrium described in
% "Explosive Instability and Erupting
% Flux Tubes in a Magnetised Plasma Atmosphere" by Cowley, Cowley, 
% Henneberg, Wilson, arXiv:1411.7797v1,  Eq. (54) 
% see Eq. (57) and (58)
% input parameters are: x_rho, x_B, A, B1=B1^2, B2=B2^2, q, x0




function [gamma, C4] = DriveCoefficients(A, B1, B2, x0, xB, xrho)

    gamma = -(A^2 * B1 - (A^2 * B2)./(cosh(x0-xB).^2))* pi^2 ...
        - sinh(x0-xrho)./(cosh(x0-xrho).^3);
    C4 = 8/(3*pi)*(3*tanh(x0-xrho).^2-1.0)./(cosh(x0-xrho).^2);
end

    