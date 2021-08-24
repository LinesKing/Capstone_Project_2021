function [] = iprint(bp,iter,x,lambda,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU,auxdata)
%%%% This function is print the iteration
%%%% Input: bp, iter, x, lam, zL, zU, alphaPr, alphaDu, s, xL, xU, bL, bU
%%%%    bp - barrier problem
%%%%    iter - iteration number
%%%%    x - primal variable
%%%%    lambda - lambda (lagragian multipliers)
%%%%    zL - lower bound for z
%%%%    zU - upper bound for z
%%%%    alphaPr - alpha in line search for primal variables
%%%%    alphaDu - alpha in line search for primal variables
%%%%    s - slack variable
%%%%    xL: lower bound for x
%%%%    xU: upper bound for x
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: Na

    screen = 1;
    if (mod(iter,15)==0)
       fprintf(screen,'   Iter       Merit   Objective   log10(mu)        Pcvg        Dcvg    alphaPr    alphaDu\n');
    end
    du = sum(abs(objectiveGradient(x,s,auxdata)' + constraintJacobian(x,bL,bU,auxdata)'*lambda' - zL' + zU'));
    me = merit(bp,x,xL,xU,s,bL,bU,auxdata);
    ob = objective(x,auxdata);
    fprintf(screen,'  %5i %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n',iter,me,ob,log10(bp.mu),residualAbsSum(x,s,bL,bU,auxdata),du,alphaPr,alphaDu);
end