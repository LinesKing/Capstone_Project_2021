function [J] = constraintJacobian(chi,bL,bU,auxdata)
%%%% This function is the rearrange Jacobian matrix
%%%% Input: x, bL, bU
%%%%    x - primal variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%%    auxdata - the auxiliary data
%%%% Output: pd
%%%%    J: jacobian matrix

    % Define Constants
    [N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner, A, b] = deal(auxdata{:});

    % Number of constraints
    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraintsObject = 2;
    nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

    % Allocate space for Jacobian Sparse
    J = spalloc(nConstraints*N,nStates*(N+1),N*32); 

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
% % %     J = double(subs(jacobian(c, xVector),xVector,chi'));
    
    % Counters
    k = 1; % Column
    m = 1; % Row

    % Fill in Jacobian Matrix
    for i = 1:N
        % Jacobian of System Dynamics
        J(m,k) = -1; J(m,k+2) = -chi(k+10)*cos(chi(k+3)); J(m,k+3) = chi(k+10)*chi(k+2)*sin(chi(k+3)); J(m,k+nStates) = 1; J(m,k+10) = -chi(k+2)*cos(chi(k+3));
        J(m+1,k+1) = -1; J(m+1,k+2) = -chi(k+10)*sin(chi(k+3)); J(m+1,k+3) = -chi(k+10)*chi(k+2)*cos(chi(k+3)); J(m+1,k+1+nStates) = 1; J(m+1,k+10) = -chi(k+2)*sin(chi(k+3)); 
        J(m+2,k+2) = -1; J(m+2,k+4) = -chi(k+10); J(m+2,k+2+nStates) = 1; J(m+2,k+10) = -chi(k+4);
        J(m+3,k+3) = -1; J(m+3,k+5) = -chi(k+10); J(m+3,k+3+nStates) = 1; J(m+3,k+10) = -chi(k+5);
        m = m + nConstraintsSystem;

        % Jacobian of Sticking to Triplet Constraint
        if (xInner(i) ~= xOuter(i))
            J(m,k) = -(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i));
            J(m,k+1) = 1;
        end
        m = m + nConstraintsTriplet;

        % Jacobian of Obstacle Constraints
        J(m, k) = chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1);
        J(m, k+1) = chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2);
        J(m, k+6) = A(1,1)*chi(k) + A(1,2)*chi(k+1) - b(1,1);
        J(m, k+7) = A(2,1)*chi(k) + A(2,2)*chi(k+1) - b(2,1);
        J(m, k+8) = A(3,1)*chi(k) + A(3,2)*chi(k+1) - b(3,1);
        J(m, k+9) = A(4,1)*chi(k) + A(4,2)*chi(k+1) - b(4,1);

        J(m+1, k+6) = 2*A(1,1)*(chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1)) + 2*A(1,2)*(chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2));
        J(m+1, k+7) = 2*A(2,1)*(chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1)) + 2*A(2,2)*(chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2));
        J(m+1, k+8) = 2*A(3,1)*(chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1)) + 2*A(3,2)*(chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2));
        J(m+1, k+9) = 2*A(4,1)*(chi(k+6)*A(1,1) + chi(k+7)*A(2,1) + chi(k+8)*A(3,1) + chi(k+9)*A(4,1)) + 2*A(4,2)*(chi(k+6)*A(1,2) + chi(k+7)*A(2,2) + chi(k+8)*A(3,2) + chi(k+9)*A(4,2));

        m = m + nConstraintsObject;
    
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