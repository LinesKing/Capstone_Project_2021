function [h] = hes(x,s,lam)
%%%% This function is to calculate the Hessian matrix (2nd derivatives)
%%%% Input: x, s, lam
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    lam: lambda (lagragian multipliers)
%%%% Output: h
%%%%    h: hessian matrix

    n = size(x,2);
    ns = size(s,2);

    h = zeros(2,2);
    h(1,1) = 2 * lam(2);
    h(1,2) = 1 + lam(1);
    h(2,1) = 1 + lam(1);
    h(2,2) = 2 * lam(2);        


    % expand hessian for inequality slacks
    for i = 1:ns
        h(n+i,n+i) = 0.0;
    end
end