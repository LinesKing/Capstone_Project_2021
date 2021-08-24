function [me] = merit(bp,x,xL,xU,s,bL,bU,auxdata)
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

    phi = barrierObjective(bp,x,xL,xU,s,bL,bU,auxdata);
    r = constaintResidual(x,s,bL,bU,auxdata);
    me = phi + bp.nu*sum(abs(r));
end