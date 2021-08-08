function [ph] = phi(bp,x,xL,xU,s,bL,bU)
%%%% This function is to calculate the lagragian function
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

    ph = obj(x);
    n = size(x,2);
    m = size(bL,2);
    for i = 1:n
       ph = ph - bp.mu * (log(x(i)-xL(i)) + log(xU(i)-x(i)));
    end
    j = 0;
    for i = 1:m
       if(bU(i)>bL(i))
          j = j + 1;
          ph = ph - bp.mu * (log(s(j)-bL(i)) + log(bU(i)-s(j)));
       end
    end
end