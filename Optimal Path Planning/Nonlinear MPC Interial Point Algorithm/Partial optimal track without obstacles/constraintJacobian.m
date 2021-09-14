function [J] = constraintJacobian(x,bL,bU,auxdata)
%%%% This function is the rearrange Jacobian matrix
%%%% Input: x, bL, bU
%%%%    x - primal variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%%    auxdata - the auxiliary data
%%%% Output: pd
%%%%    J: jacobian matrix

    % Define Constants
    [N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner] = deal(auxdata{:});

    % Number of constraints
    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraints = nConstraintsTriplet + nConstraintsSystem;

    % Allocate space for Jacobian Sparse
    J = spalloc(nConstraints*N,nStates*(N+1),N*18); 

    % Counters
    k = 1; % Column
    m = 1; % Row

    % Fill in Jacobian Matrix
    for i = 1:N
        % Jacobian of System Dynamics
        J(m,k) = -1; J(m,k+2) = -x(k+6); J(m,k+nStates) = 1; J(m,k+6) = -x(k+2);
        J(m+1,k+1) = -1; J(m+1,k+3) = -x(k+6); J(m+1,k+1+nStates) = 1; J(m+1,k+6) = -x(k+3);
        J(m+2,k+2) = -1; J(m+2,k+4) = -x(k+6); J(m+2,k+2+nStates) = 1; J(m+2,k+6) = -x(k+4);
        J(m+3,k+3) = -1; J(m+3,k+5) = -x(k+6); J(m+3,k+3+nStates) = 1; J(m+3,k+6) = -x(k+5);
        m = m + nConstraintsSystem;

        % Jacobian of Sticking to Triplet Constraint
        if (xInner(i) ~= xOuter(i))
            J(m,k) = -(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i));
            J(m,k+1) = 1;
        end
        m = m + nConstraintsTriplet;

        % Next Iteration
        k = k + nStates;
    end

    m = size(J,1);
    n = size(J,2);

    % add slack variables for inequality constraints
    k = 0;
    for i = 1:m
       if(bU(i)>bL(i))
           k = k + 1;
           J(i,n+k) = -1;
       end
    end
end