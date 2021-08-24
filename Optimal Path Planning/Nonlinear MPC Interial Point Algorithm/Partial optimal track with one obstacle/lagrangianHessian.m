function [h] = lagrangianHessian(x,s,lam,auxdata)
%%%% This function is to calculate the Hessian matrix (2nd derivatives)
%%%% Input: x, s, lam
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    lam: lambda (lagragian multipliers)
%%%% Output: h
%%%%    h: hessian matrix

    % Define constants
    [N, nStates, weightDist, weightVel, weightTheta, ~, ~, ~, ~, ~, A, ~] = deal(auxdata{:});
    
    n = size(x,2);
    ns = size(s,2);
    
    % Number of constraints
    nConstraintsSystem = 4;
    nConstraintsTriplet = 1;
    nConstraintsObject = 2;
    nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

    % Generate gradadient of objective function:  
    h = zeros((N+1)*nStates,(N+1)*nStates);

    % Counter
    k = 1; % lam
    m = 1; % x
    
    % Get upper triangular part of h
    for i = 1:N+1
        if (i == 1)
            h(m,m) = 2*weightDist; h(m,m+nStates) = -2*weightDist; h(m,m+6) = lam(k+5)*A(1,1); h(m,m+7) = lam(k+5)*A(2,1); h(m,m+8) = lam(k+5)*A(3,1); h(m,m+9) = lam(k+5)*A(4,1); 
            h(m+1,m+1) = 2*weightDist; h(m+1,m+1+nStates) = -2*weightDist; h(m,m+6) = lam(k+5)*A(1,2); h(m,m+7) = lam(k+5)*A(2,2); h(m,m+8) = lam(k+5)*A(3,2); h(m,m+9) = lam(k+5)*A(4,2);
            h(m+2,m+2) = 2*weightVel; h(m+2,m+2+nStates) = -2*weightVel; h(m+2,m+3) = lam(k)*x(m+10)*sin(x(m+3)) - lam(k+1)*x(m+10)*cos(x(m+3)); h(m+2,m+10) = -lam(k)*cos(x(m+3))-lam(k+1)*sin(x(m+3));
            h(m+3,m+3) = 2*weightTheta + lam(k)*x(m+10)*x(m+2)*cos(x(m+3)) + lam(k+1)*x(m+10)*x(m+2)*sin(x(m+3)); h(m+3,m+3+nStates) = -2*weightVel; h(m+3,m+10) = lam(k)*x(m+2)*sin(x(m+3))-lam(k+1)*x(m+2)*cos(x(m+3));
            h(m+4,m+10) = -lam(k+2);
            h(m+5,m+10) = -lam(k+3);
            h(m+6,m+6) = 2*lam(k+6)*(A(1,1)^2+A(1,2)^2); h(m+6,m+7) = 2*lam(k+6)*(A(1,1)*A(2,1)+A(1,2)*A(2,2)); h(m+6,m+8) = 2*lam(k+6)*(A(1,1)*A(3,1)+A(1,2)*A(3,2)); h(m+6,m+9) = 2*lam(k+6)*(A(1,1)*A(4,1)+A(1,2)*A(4,2));
            h(m+7,m+7) = 2*lam(k+6)*(A(2,1)^2+A(2,2)^2); h(m+7,m+8) = 2*lam(k+6)*(A(2,1)*A(3,1)+A(2,2)*A(3,2)); h(m+7,m+9) = 2*lam(k+6)*(A(2,1)*A(4,1)+A(2,2)*A(4,2));
            h(m+8,m+8) = 2*lam(k+6)*(A(3,1)^2+A(3,2)^2); h(m+8,m+9) = 2*lam(k+6)*(A(3,1)*A(4,1)+A(3,2)*A(4,2));
            h(m+9,m+9) = 2*lam(k+6)*(A(4,1)^2+A(4,2)^2);
        elseif (i == N+1)
            h(m,m) = 2*weightDist;
            h(m+1,m+1) = 2*weightDist;
            h(m+2,m+2) = 2*weightVel;
            h(m+3,m+3) = 2*weightTheta;
        else
            h(m,m) = 4*weightDist; h(m,m+nStates) = -2*weightDist; h(m,m+6) = lam(k+5)*A(1,1); h(m,m+7) = lam(k+5)*A(2,1); h(m,m+8) = lam(k+5)*A(3,1); h(m,m+9) = lam(k+5)*A(4,1); 
            h(m+1,m+1) = 4*weightDist; h(m+1,m+1+nStates) = -2*weightDist; h(m,m+6) = lam(k+5)*A(1,2); h(m,m+7) = lam(k+5)*A(2,2); h(m,m+8) = lam(k+5)*A(3,2); h(m,m+9) = lam(k+5)*A(4,2);
            h(m+2,m+2) = 4*weightVel; h(m+2,m+2+nStates) = -2*weightVel; h(m+2,m+3) = lam(k)*x(m+10)*sin(x(m+3)) - lam(k+1)*x(m+10)*cos(x(m+3)); h(m+2,m+10) = -lam(k)*cos(x(m+3))-lam(k+1)*sin(x(m+3));
            h(m+3,m+3) = 4*weightTheta + lam(k)*x(m+10)*x(m+2)*cos(x(m+3)) + lam(k+1)*x(m+10)*x(m+2)*sin(x(m+3)); h(m+3,m+3+nStates) = -2*weightVel; h(m+3,m+10) = lam(k)*x(m+2)*sin(x(m+3))-lam(k+1)*x(m+2)*cos(x(m+3));
            h(m+4,m+10) = -lam(k+2);
            h(m+5,m+10) = -lam(k+3);
            h(m+6,m+6) = 2*lam(k+6)*(A(1,1)^2+A(1,2)^2); h(m+6,m+7) = 2*lam(k+6)*(A(1,1)*A(2,1)+A(1,2)*A(2,2)); h(m+6,m+8) = 2*lam(k+6)*(A(1,1)*A(3,1)+A(1,2)*A(3,2)); h(m+6,m+9) = 2*lam(k+6)*(A(1,1)*A(4,1)+A(1,2)*A(4,2));
            h(m+7,m+7) = 2*lam(k+6)*(A(2,1)^2+A(2,2)^2); h(m+7,m+8) = 2*lam(k+6)*(A(2,1)*A(3,1)+A(2,2)*A(3,2)); h(m+7,m+9) = 2*lam(k+6)*(A(2,1)*A(4,1)+A(2,2)*A(4,2));
            h(m+8,m+8) = 2*lam(k+6)*(A(3,1)^2+A(3,2)^2); h(m+8,m+9) = 2*lam(k+6)*(A(3,1)*A(4,1)+A(3,2)*A(4,2));
            h(m+9,m+9) = 2*lam(k+6)*(A(4,1)^2+A(4,2)^2);
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