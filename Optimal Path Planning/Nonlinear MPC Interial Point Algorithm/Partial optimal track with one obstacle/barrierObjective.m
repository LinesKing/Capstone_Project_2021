function [phi] = barrierObjective(bp,x,xL,xU,s,bL,bU,auxdata)
%%%% This function is to calculate the barrier function
%%%% Input: bp, x, xL, xU, s, bL, bU
%%%%    bp - barrier problem
%%%%    x - primal variable
%%%%    xL: lower bound for x
%%%%    xU: upper bound for x
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: ph
%%%%    ph: phi(lagragian function)

    phi = objective(x,auxdata);
    n = size(x,2);
    m = size(bL,2);
    for i = 1:n
        if (x(i)-xL(i))~=0
            phi = phi - bp.mu * (log(x(i)-xL(i)));
        end
        if (xU(i)-x(i))~=0
            phi = phi - bp.mu * (log(xU(i)-x(i)));
        end
    end
    j = 0;
    for i = 1:m
       if(bU(i)>bL(i))
          j = j + 1;
          phi = phi - bp.mu * (log(s(j)-bL(i)) + log(bU(i)-s(j)));
       end
    end
end