function [ob] = objective(chi,auxdata)
%%%% This function is to achieve the objective function value
%%%% Input: x
%%%%    x - primal variable
%%%% Output: ob
%%%%    ob: object

    % Define constants
    [N, nStates, weightDist, weightVel, weightTheta, weightTime] = deal(auxdata{1:6});

    % Generate objective function:
    %   \Sigma_{i=1}^N (\delta_x_i^T W_d \delta_x_i + 
    %   \delta_y_i^T W_d \delta_y_i + \delta_v_x_i^T W_v 
    %   \delta_v_x_i + \delta_v_y_i^T W_v \delta_v_y_i + W_t t_i)
    
% % %     % Initialize the variables
% % %     x = sym('x', [N+1 nStates]);
% % %     dx = x(2:40,:) - x(1:39,:);
% % %     
% % %     f = weightDist*dx(:,1)'*dx(:,1) + weightDist*dx(:,2)'*dx(:,2)...
% % %         + weightVel*dx(:,3)'*dx(:,3) + weightTheta*dx(:,4)'*dx(:,4)...
% % %         + weightTime*sum(x(:,11));
% % %     
% % %     xVector = x(:);
% % %     ob = double(subs(f,xVector,chi'));
    
    ob = 0;
    k = 1;
    for i = 1:N
        ob = ob + weightDist*(chi(k) - chi(k+nStates))^2 ...
                + weightDist*(chi(k+1) - chi(k+1+nStates))^2 ...
                + weightVel*(chi(k+2) - chi(k+2+nStates))^2 ...
                + weightTheta*(chi(k+3) - chi(k+3+nStates))^2 ...
                + weightTime*chi(k+10);

        % Next interation
        k = k + nStates;
    end
end