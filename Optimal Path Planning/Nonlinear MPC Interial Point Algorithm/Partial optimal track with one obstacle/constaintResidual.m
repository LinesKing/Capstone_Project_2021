function [r] = constaintResidual(chi,s,bL,bU,auxdata)
%%%% This function is residual with equality constants and inequality slack variables
%%%% Input: x, s, bL, bU
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: r
%%%%    r: residual

    % Define constants
    [N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner, A, b] = deal(auxdata{:});

    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraintsObject = 2;
    nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

% % %     % Initialize the variables
% % %     x = sym('x', [N+1 nStates]);
% % %     c = sym('c', [N*nConstraints 1]);
% % %     m = 1; % Row
% % %     
% % %     for i = 1:N
% % %         % System dynamics constraints
% % %         c(m, 1) = x(i+1,1) - x(i,1) - x(i,11)*x(i,3)*cos(x(i,4));
% % %         c(m+1, 1) = x(i+1,2) - x(i,2) - x(i,11)*x(i,3)*sin(x(i,4));
% % %         c(m+2, 1) = x(i+1,3) - x(i,3) - x(i,11)*x(i,5);
% % %         c(m+3, 1) = x(i+1,4) - x(i,4) - x(i,11)*x(i,6);
% % %         m = m + nConstraintsSystem;
% % % 
% % %         % Constraints stick to triplet (linear interpolation)
% % %         if (xInner(i) == xOuter(i) || i == 1)
% % %             c(m,1) = 0;
% % %         else
% % %             c(m,1) = x(i,2) - x(i,1)*(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i)) - (yInner(i)-(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i))*xInner(i)); 
% % %         end
% % %         m = m + nConstraintsTriplet;
% % %     
% % %         % Obstacle Constraints
% % %         %   (A*p_i - b)^T \lambda > d_{min}
% % %         p = [x(i,1) x(i,2)]';
% % %         lambda = [x(i,7) x(i,8) x(i,9) x(i,10)]';
% % %         c(m,1) = (A*p-b)'*lambda;
% % %         %   ||A^T \lambda|| <= 1
% % %         c(m+1,1) = (A(:,1)'*lambda)^2 + (A(:,2)'*lambda)^2; 
% % %         m = m + nConstraintsObject;      
% % %     end
% % % 
% % %     xVector = x(:);
% % %     r = double(subs(c,xVector,chi'));
    
    % Counter
    k = 1; % Column
    m = 1; % Row
    for i = 1:N
        % System dynamics constraints
        r(m, 1) = chi(k+nStates) - chi(k) - chi(k+10)*chi(k+2)*cos(chi(k+3));
        r(m+1, 1) = chi(k+1+nStates) - chi(k+1) - chi(k+10)*chi(k+2)*sin(chi(k+3));
        r(m+2, 1) = chi(k+2+nStates) - chi(k+2) - chi(k+10)*chi(k+4);
        r(m+3, 1) = chi(k+3+nStates) - chi(k+3) - chi(k+10)*chi(k+5);
        m = m + nConstraintsSystem;

        % Constraints stick to triplet (linear interpolation)
        if (xInner(i) == xOuter(i) || i == 1)
            r(m,1) = 0;
        else
            r(m,1) = chi(k+1) - chi(k)*(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i)) - (yInner(i)-(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i))*xInner(i)); 
        end
        m = m + nConstraintsTriplet;
    
        % Obstacle Constraints
        %   (A*p_i - b)^T \lambda > d_{min}
        r(m,1) = chi(k+6)*(A(1,1)*chi(k) + A(1,2)*chi(k+1) - b(1,1))+ chi(k+7)*(A(2,1)*chi(k) + A(2,2)*chi(k+1) - b(2,1)) + chi(k+8)*(A(3,1)*chi(k) + A(3,2)*chi(k+1) - b(3,1)) + chi(k+9)*(A(4,1)*chi(k) + A(4,2)*chi(k+1) - b(4,1));
        %   ||A^T \lambda|| <= 1
        r(m+1,1) = (chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1))^2 + (chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2))^2; 
        m = m + nConstraintsObject;   
        
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

