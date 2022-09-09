% Fourth Order Runge-Kutta solver (f,x0,y0,h,N)
% e.g. see Advanced Engineering Mathematics (10th edition)
% page 902
% inputs(f, initial location x0, initial values of y' = y0, step size h,
% number of steps N)
% with y'=f(x,y) and y(x0)=y0

function [x,y]= RungeKuttaSolver(f,x0, y0, h, N)
    x(1)=x0;
    y(1)=y0;
    for i = 1:N
        k1=h*f(x(i),y(i));
        k2=h*f(x(i)+0.5*h, y(i)+0.5*k1);
        k3=h*f(x(i)+0.5*h, y(i)+0.5*k2);
        k4=h*f(x(i)+h,y(i)+k3);
        
        x(i+1)=x(i)+h;
        y(i+1)=y(i)+1/6*(k1+2*k2+2*k3+k4);
    end
    
end