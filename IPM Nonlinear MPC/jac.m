function [pd] = jac(x,bL,bU)
%%%% This function is the rearrange Jacobian matrix
%%%% Input: x, bL, bU
%%%%    x - primal variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: pd
%%%%    pd: jacobian matrix

    pd(1,1) = x(2);
    pd(1,2) = x(1);
    pd(2,1) = 2*x(1);
    pd(2,2) = 2*x(2);

    m = size(pd,1);
    n = size(pd,2);

    % add slack variables for inequality constraints
    k = 0;
    for i = 1:m
       if(bU(i)>bL(i))
           k = k + 1;
           pd(i,n+k) = -1;
       end
    end
end