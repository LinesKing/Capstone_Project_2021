function [h] = lagrangianHessian(x,s,lam,auxdata)
%%%% This function is to calculate the Hessian matrix (2nd derivatives)
%%%% Input: x, s, lam
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    lam: lambda (lagragian multipliers)
%%%% Output: h
%%%%    h: hessian matrix

    % Define constants
    [N, nStates, weightDist, weightVel] = deal(auxdata{1:4});
    
    n = size(x,2);
    ns = size(s,2);
    
    % Number of constraints
    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraints = nConstraintsTriplet + nConstraintsSystem;

    % Generate gradadient of objective function:  
    h = zeros((N+1)*nStates,(N+1)*nStates);

    % Counter
    k = 1; % lam
    m = 1; % x
    
    % Get upper triangular part of h
    for i = 1:N+1
        if (i == 1)
            h(m,m) = 2*weightDist; h(m,m+nStates) = -2*weightDist;
            h(m+1,m+1) = 2*weightDist; h(m+1,m+1+nStates) = -2*weightDist;
            h(m+2,m+2) = 2*weightVel; h(m+2,m+2+nStates) = -2*weightVel; h(m+2,m+6) = -lam(k);
            h(m+3,m+3) = 2*weightVel; h(m+3,m+3+nStates) = -2*weightVel; h(m+3,m+6) = -lam(k+1);
            h(m+4,m+6) = -lam(k+2);
            h(m+5,m+6) = -lam(k+3);
        elseif (i == N+1)
            h(m,m) = 2*weightDist;
            h(m+1,m+1) = 2*weightDist;
            h(m+2,m+2) = 2*weightVel;
            h(m+3,m+3) = 2*weightVel;
        else
            h(m,m) = 4*weightDist; h(m,m+nStates) = -2*weightDist;
            h(m+1,m+1) = 4*weightDist; h(m+1,m+1+nStates) = -2*weightDist;
            h(m+2,m+2) = 4*weightVel; h(m+2,m+2+nStates) = -2*weightVel; h(m+2,m+6) = -lam(k);
            h(m+3,m+3) = 4*weightVel; h(m+3,m+3+nStates) = -2*weightVel; h(m+3,m+6) = -lam(k+1);
            h(m+4,m+6) = -lam(k+2);
            h(m+5,m+6) = -lam(k+3);
        end
        
        m = m + nStates;
        k = k + nConstraints;
        
        ns = size(s,2);
    end    

    % Get full h
    h = h' + triu(h,1);
    
    % expand hessian for inequality slacks
    for i = 1:ns
        h(n+i,n+i) = 0.0;
    end
end