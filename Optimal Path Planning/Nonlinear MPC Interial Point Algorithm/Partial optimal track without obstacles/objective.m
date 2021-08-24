function [ob] = objective(x,auxdata)
%%%% This function is to achieve the objective function value
%%%% Input: x
%%%%    x - primal variable
%%%% Output: ob
%%%%    ob: object

    % Define constants
    [N, nStates, weightDist, weightVel, weightTime] = deal(auxdata{1:5});

    % Generate objective function:
    %   \Sigma_{i=1}^N (\delta_x_i^T W_d \delta_x_i + 
    %   \delta_y_i^T W_d \delta_y_i + \delta_v_x_i^T W_v 
    %   \delta_v_x_i + \delta_v_y_i^T W_v \delta_v_y_i + W_t t_i)
    ob = 0;
    k = 1;
    for i = 1:N
        ob = ob + weightDist*(x(k) - x(k+nStates))^2 ...
              + weightDist*(x(k+1) - x(k+1+nStates))^2 ...
              + weightVel*(x(k+2) - x(k+2+nStates))^2 ...
              + weightVel*(x(k+3) - x(k+3+nStates))^2 ...
              + weightTime*x(k+6);

        % Next interation
        k = k + nStates;
    end
end