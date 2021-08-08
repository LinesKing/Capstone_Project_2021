function [me] = merit(bp,x,xL,xU,s,bL,bU)
%%%% This function is the merit function
%%%% Input: bp, x, xL, xU, s, bL, bU
%%%%    bp - barrier problem
%%%%    x - primal variable
%%%%    xL: lower bound for x
%%%%    xU: upper bound for x
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: me
%%%%    me: merit

    ph = phi(bp,x,xL,xU,s,bL,bU);
    r = res(x,s,bL,bU);
    me = ph + bp.nu*sum(abs(r));
end