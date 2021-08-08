function [] = iprint(bp,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU)
screen = 1;
if (mod(iter,15)==0)
   fprintf(screen,'   Iter       Merit   Objective   log10(mu)        Pcvg        Dcvg    alphaPr    alphaDu\n');
end
du = sum(abs(objGrad(x,s)' + jac(x,bL,bU)'*lam' - zL' + zU'));
me = merit(bp,x,xL,xU,s,bL,bU);
ob = obj(x);
fprintf(screen,'  %5i %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n',iter,me,ob,log10(bp.mu),theta(x,s,bL,bU),du,alphaPr,alphaDu);
