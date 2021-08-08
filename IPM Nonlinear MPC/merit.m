function [me] = merit(bp,x,xL,xU,s,bL,bU)
    ph = phi(bp,x,xL,xU,s,bL,bU);
    r = res(x,s,bL,bU);
    me = ph + bp.nu*sum(abs(r));