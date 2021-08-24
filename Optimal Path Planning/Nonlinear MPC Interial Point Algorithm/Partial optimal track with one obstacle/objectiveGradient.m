function [g] = objectiveGradient(chi,s,auxdata)
%%%% This function is to obtain objective function gradient
%%%% Input: x, s
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    auxdata - the auxiliary data
%%%% Output: g
%%%%    g: gradient of object function

    % Define constants
    [N, nStates, weightDist, weightVel, weightTheta, weightTime] = deal(auxdata{1:6});

% % %     % Initialize the variables
% % %     x = sym('x', [N+1 nStates]);
% % %     dx = x(2:40,:) - x(1:39,:);
% % %     
% % %     f = weightDist*dx(:,1)'*dx(:,1) + weightDist*dx(:,2)'*dx(:,2)...
% % %         + weightVel*dx(:,3)'*dx(:,3) + weightTheta*dx(:,4)'*dx(:,4)...
% % %         + weightTime*sum(x(:,11));
% % %     
% % %     xVector = x(:);
% % %     g = double(subs(gradient(f,xVector),xVector,chi'));
    
    
    % Generate gradadient of objective function:  
    g = zeros((N+1)*nStates,1);

    
    % Counter
    k = 1; % Column
    m = 1; % Row

    for i = 1:N+1
        if (i == 1)
            g(m,1) = 2*weightDist*(chi(k)-chi(k+nStates));
            g(m+1,1) = 2*weightDist*(chi(k+1)-chi(k+1+nStates));
            g(m+2,1) = 2*weightVel*(chi(k+2)-chi(k+2+nStates));
            g(m+3,1) = 2*weightTheta*(chi(k+3)-chi(k+3+nStates));
            g(m+10,1) = weightTime;
        elseif (i == N+1)
            g(m,1) = -2*weightDist*(chi(k-nStates)-chi(k));
            g(m+1,1) = -2*weightDist*(chi(k+1-nStates)-chi(k+1));
            g(m+2,1) = -2*weightVel*(chi(k+2-nStates)-chi(k+2));
            g(m+3,1) = -2*weightTheta*(chi(k+3-nStates)-chi(k+3));
            g(m+10,1) = weightTime;
        else
            g(m,1) = -2*weightDist*(chi(k-nStates)-chi(k)) + 2*weightDist*(chi(k)-chi(k+nStates));
            g(m+1,1) = -2*weightDist*(chi(k+1-nStates)-chi(k+1)) + 2*weightDist*(chi(k+1)-chi(k+1+nStates));
            g(m+2,1) = -2*weightVel*(chi(k+2-nStates)-chi(k+2)) + 2*weightVel*(chi(k+2)-chi(k+2+nStates));
            g(m+3,1) = -2*weightTheta*(chi(k+3-nStates)-chi(k+3)) + 2*weightTheta*(chi(k+3)-chi(k+3+nStates));
            g(m+10,1) = weightTime;
        end
        m = m + nStates;

        % Next Interation
        k = k + nStates;
        
        ns = size(s,2);
    end

    ns = size(s,2);
    g = [g' zeros(1,ns)];
end