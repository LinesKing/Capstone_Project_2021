function [ob] = objective(x,auxdata)
%%%% This function is to achieve the objective function value
%%%% Input: x
%%%%    x - primal variable
%%%% Output: ob
%%%%    ob: object

    % Define constants
    [N, nStates, weightDist, weightVel, weightAcc, weightTime] = deal(auxdata{1:6});

    % Generate objective function:
    ob = 0;
    k = 1;
    for i = 1:N
        ob = ob + weightDist*(x(k) - x(k+nStates))^2 ...
              + weightDist*(x(k+1) - x(k+1+nStates))^2 ...
              + weightVel*(x(k+2) - x(k+2+nStates))^2 ...
              + weightVel*(x(k+3) - x(k+3+nStates))^2 ...
              + weightAcc*(x(k+4) - x(k+4+nStates))^2 ...
              + weightAcc*(x(k+5) - x(k+5+nStates))^2 ...
              + weightTime*x(k+6);

        % Next interation
        k = k + nStates;
    end
end