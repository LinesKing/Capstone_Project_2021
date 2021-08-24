function [g] = objectiveGradient(x,s,auxdata)
%%%% This function is to obtain objective function gradient
%%%% Input: x, s
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    auxdata - the auxiliary data
%%%% Output: g
%%%%    g: gradient of object function

    % Define constants
    [N, nStates, weightDist, weightVel, weightTime] = deal(auxdata{1:5});


    % Generate gradadient of objective function:  
    g = zeros((N+1)*nStates,1);

    % Counter
    k = 1; % Column
    m = 1; % Row

    for i = 1:N+1
        if (i == 1)
            g(m,1) = 2*weightDist*(x(k)-x(k+nStates));
            g(m+1,1) = 2*weightDist*(x(k+1)-x(k+1+nStates));
            g(m+2,1) = 2*weightVel*(x(k+2)-x(k+2+nStates));
            g(m+3,1) = 2*weightVel*(x(k+3)-x(k+3+nStates));
            g(m+6,1) = weightTime;
        elseif (i == N+1)
            g(m,1) = -2*weightDist*(x(k-nStates)-x(k));
            g(m+1,1) = -2*weightDist*(x(k+1-nStates)-x(k+1));
            g(m+2,1) = -2*weightVel*(x(k+2-nStates)-x(k+2));
            g(m+3,1) = -2*weightVel*(x(k+3-nStates)-x(k+3));
            g(m+6,1) = weightTime;
        else
            g(m,1) = -2*weightDist*(x(k-nStates)-x(k)) + 2*weightDist*(x(k)-x(k+nStates));
            g(m+1,1) = -2*weightDist*(x(k+1-nStates)-x(k+1)) + 2*weightDist*(x(k+1)-x(k+1+nStates));
            g(m+2,1) = -2*weightVel*(x(k+2-nStates)-x(k+2)) + 2*weightVel*(x(k+2)-x(k+2+nStates));
            g(m+3,1) = -2*weightVel*(x(k+3-nStates)-x(k+3)) + 2*weightVel*(x(k+3)-x(k+3+nStates));
            g(m+6,1) = weightTime;
        end
        m = m + nStates;

        % Next Interation
        k = k + nStates;
        
        ns = size(s,2);
    end
    g = [g' zeros(1,ns)];
end