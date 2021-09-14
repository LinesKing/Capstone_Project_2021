function [r] = constaintResidual(x,s,bL,bU,auxdata)
%%%% This function is residual with equality constants and inequality slack variables
%%%% Input: x, s, bL, bU
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: r
%%%%    r: residual

    % Define constants
    [N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner] = deal(auxdata{:});

    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraints = nConstraintsTriplet + nConstraintsSystem;

    % Generate matrix of constraint functions
    r = zeros(nConstraints*N,1);

    % Counter
    k = 1; % Column
    m = 1; % Row
    for i = 1:N
        % System dynamics constraints
        r(m, 1) = x(k+nStates) - x(k) - x(k+6)*x(k+2);
        r(m+1, 1) = x(k+1+nStates) - x(k+1) - x(k+6)*x(k+3);
        r(m+2, 1) = x(k+2+nStates) - x(k+2) - x(k+6)*x(k+4);
        r(m+3, 1) = x(k+3+nStates) - x(k+3) - x(k+6)*x(k+5);
        m = m + nConstraintsSystem;

        % Constraints stick to triplet (linear interpolation)
        if (xInner(i) == xOuter(i) || i == 1)
            r(m,1) = 0;
        else
            r(m,1) = x(k+1) - x(k)*(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i)) - (yInner(i)-(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i))*xInner(i)); 
        end
        m = m + nConstraintsTriplet;

        % Next Interation
        k = k + nStates;
    end
    
    j = 0;
    for i = 1:size(r,2)
        if (bU(i)==bL(i))
            % equality constant
            r(i) = r(i) - bL(i);
        else
            % inequality slack
            j = j + 1;
            r(i) = r(i) - s(j);
        end
    end
end

