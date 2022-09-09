% SecantMethod: root finding algorithm

function[xroot]= SecantMethod(f,x0,x1,itr)
    success=0;
    for i=1:itr
        x2=(x0*f(x1)-x1*f(x0))/(f(x1)-f(x0));
        if abs(x2-x1)<tol
            xroot=x2;
            success=1;
            break;
        else
            x0=x1;
            x1=x2;
        end
    end
    xroot=x2;
end