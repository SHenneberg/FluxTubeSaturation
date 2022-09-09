% Runge-Kutta-Nystroem solver (f,x0,y0, dy0,h,N)
% e.g. see Advanced Engineering Mathematics (10th edition)
% page 916
% inputs(f, initial location x0, initial values of y = y0, and y'=dy0, step size h,
% number of steps N)
% with y''=f(x,y,dy) and y(x0)=y0, dy(x0)=dy0

function [x,y,dy]= RungeKuttaNystroemSolver(f,x0, y0,dy0, h, N)
    x(1)=x0;
    y(1)=y0;
    dy(1)=dy0;
    for i = 1:N
        k1=0.5*h*f(x(i),y(i),dy(i));
        K=0.5*h*(dy(i)+0.5*k1);
        k2=0.5*h*f(x(i)+0.5*h, y(i)+K, dy(i)+k1);
        k3=0.5*h*f(x(i)+0.5*h, y(i)+K, dy(i)+k2);
        L=h*(dy(i)+k3);
        k4=0.5*h*f(x(i)+h, y(i)+L, dy(i)+2*k3);
        
        x(i+1)=x(i)+h;
        y(i+1)=y(i)+h*(dy(i)+1/3*(k1+k2+k3));
        dy(i+1)=dy(i)+1/3*(k1+2*k2+2*k3+k4);
    end
    
end