function [x] = OBCAOneObsTimeOptimalPathPlanning(xInner, yInner, xOuter, yOuter, warm, initPoint, finalPoint, object, dMin)
%%%% This function is to find the optimal path with OBCA method.
%%%% Input:
%%%% 	xInner - Inner bound of x
%%%%    yInner - Inner bound of y
%%%%    xOuter - Outer bound of x
%%%%    yOuter - Outer bound of y
%%%%    warm - warm starting points
%%%%    initPoint - Initial points
%%%% Output:
%%%%    sol - Optimal solution

% Receeding horizon
N = length(xOuter) - 1;

% Max time between points (Must be chosen to something physcally poissble)
tMax = 0.05; 

% Define state/input constraints
vMin = 0;
thetaMin = -inf;
aMin = -5;
omegaMin = -2*pi;
vMax = 2;
thetaMax = inf;
aMax = 5;
omegaMax = 2*pi;

% Define objetive function weightings
weightDist = 0.1;
weightVel = 0.5;
weightTheta = 0.1;
weightTime = 1;

% Define constants (number of variables)
nSystemStates = 4; % [x, y, v, theta]
nControlStates = 2; % [a, omega]
nObcaVariables = 4; % [lambda1, lambda2, lambda3, lambda4]
nTimeVariable = 1;
nStates = nSystemStates + nControlStates + nObcaVariables + nTimeVariable; % [x, y, v, theta, a, w, lambda1, lambda2, lambda3, lambda4]

% Initialise decision variables, lower bounds and upper bounds
x0 = [];
options.lb = [];
options.ub = [];
for i = 1:N+1
    % Set initial start point to center of track (warm start)
    % [x, y, v, theta, a, omega, lambda1, lambda2, lambda3, lambda4]
    x0 = [x0, warm.xWarm(i), warm.yWarm(i), zeros(1, nStates-2)];
   
    % Upper and lower bounds of decision variables (not constraint equations)
    % [x, y, v, theta, a, omega]
    if (i == 1) % satisfy initial conditions contraint
        options.lb = [options.lb, [initPoint(1) initPoint(2) vMin initPoint(3) aMin 0]]; %2D array with each row filled with column of vector current states
        options.ub = [options.ub, [initPoint(1) initPoint(2) vMax initPoint(3) aMax 0]];  
    else
        if (i == N+1) % satisfy final conditions contraint
            options.lb = [options.lb, [finalPoint(1) finalPoint(2) vMin finalPoint(3) aMin 0]]; %2D array with each row filled with column of vector current states
            options.ub = [options.ub, [finalPoint(1) finalPoint(2) vMax finalPoint(3) aMax 0]];  
        else
            yLb = min(yOuter(i),yInner(i));
            yUb = max(yOuter(i),yInner(i));
            xLb = min(xOuter(i),xInner(i));
            xUb = max(xOuter(i),xInner(i));
            options.lb = [options.lb, [xLb yLb vMin thetaMin aMin omegaMin]];
            options.ub = [options.ub, [xUb yUb vMax thetaMax aMax omegaMax]]; 
        end
    end
   
    % Obstacle Constraints (lambda > 0)
    % [lambda_1,lambda_2,lambda_3,lambda_4]
    options.lb = [options.lb, zeros(1, nObcaVariables)];
    options.ub = [options.ub, inf(1, nObcaVariables)];
    
    % Lower and Upper Bounds for sampling time
    options.lb = [options.lb, 0];
    options.ub = [options.ub, tMax];
end

% Initialise upper and lower bounds for constraints
options.cl = [];
options.cu = [];
for i = 1:N
    % System Constraints
    options.cl = [options.cl, [0 0 0 0]];
    options.cu = [options.cu, [0 0 0 0]];

    % Triplet contraints
    options.cl = [options.cl, 0];
    options.cu = [options.cu, 0];
    
    % Obstacle Constraints
    options.cl = [options.cl, [dMin -inf]];
    options.cu = [options.cu, [inf 1]];
end

% Set up the auxiliary data
options.auxdata = {N nStates weightDist weightVel weightTheta weightTime ...
                   xOuter yOuter xInner yInner object.A object.b};

% Set the IPOPT options.
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.tol         = 1e-7;
options.ipopt.max_iter = 1000; 
% options.ipopt.print_level = 0;

% The callback functions.
funcs.objective         = @objective;
funcs.constraints       = @constraints;
funcs.gradient          = @gradient;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobian_struct;

% Run IPOPT.
[sol info] = ipopt_auxdata(x0,funcs,options);
x = reshape(sol(1:(N+1)*nStates), [nStates,(N+1)])';

% ----------------------------------------------------------------------
function f = objective (x,auxdata)
% Define constants
[N, nStates, weightDist, weightVel, weightTheta, weightTime] = deal(auxdata{1:6});

% Generate objective function:
%   \Sigma_{i=1}^N (\delta_x_i^T W_d \delta_x_i + 
%   \delta_y_i^T W_d \delta_y_i + \delta_v^T W_v 
%   \delta_v + \delta_\theta^T W_{\theta} \delta_\theta + W_t t_i)
f = 0;
k = 1;
for i = 1:N
    f = f + weightDist*(x(k) - x(k+nStates))^2 ...
          + weightDist*(x(k+1) - x(k+1+nStates))^2 ...
          + weightVel*(x(k+2) - x(k+2+nStates))^2 ...
          + weightTheta*(x(k+3) - x(k+3+nStates))^2 ...
          + weightTime*x(k+10);
    
    % Next interation
    k = k + nStates;
end

% ----------------------------------------------------------------------
function g = gradient (x,auxdata)
% Define constants
[N, nStates, weightDist, weightVel, weightTheta, weightTime] = deal(auxdata{1:6});

% Generate gradadient of objective function:
%   
g = zeros((N+1)*nStates,1);

% Counter
k = 1; % Column
m = 1; % Row

for i = 1:N+1
    if (i == 1)
        g(m,1) = 2*weightDist*(x(k)-x(k+nStates));
        g(m+1,1) = 2*weightDist*(x(k+1)-x(k+1+nStates));
        g(m+2,1) = 2*weightVel*(x(k+2)-x(k+2+nStates));
        g(m+3,1) = 2*weightTheta*(x(k+3)-x(k+3+nStates));
        g(m+10,1) = weightTime;
    elseif (i == N+1)
        g(m,1) = -2*weightDist*(x(k-nStates)-x(k));
        g(m+1,1) = -2*weightDist*(x(k+1-nStates)-x(k+1));
        g(m+2,1) = -2*weightVel*(x(k+2-nStates)-x(k+2));
        g(m+3,1) = -2*weightTheta*(x(k+3-nStates)-x(k+3));
        g(m+10,1) = weightTime;
    else
        g(m,1) = -2*weightDist*(x(k-nStates)-x(k)) + 2*weightDist*(x(k)-x(k+nStates));
        g(m+1,1) = -2*weightDist*(x(k+1-nStates)-x(k+1)) + 2*weightDist*(x(k+1)-x(k+1+nStates));
        g(m+2,1) = -2*weightVel*(x(k+2-nStates)-x(k+2)) + 2*weightVel*(x(k+2)-x(k+2+nStates));
        g(m+3,1) = -2*weightTheta*(x(k+3-nStates)-x(k+3)) + 2*weightTheta*(x(k+3)-x(k+3+nStates));
        g(m+10,1) = weightTime;
    end
    m = m + nStates;
    
    % Next Interation
    k = k + nStates;
end

% ----------------------------------------------------------------------
function c = constraints (x,auxdata)
% Define constants
[N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner, A, b] = deal(auxdata{:});

nConstraintsSystem = 4;
nConstraintsTriplet = 1;
nConstraintsObject = 2;
nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

% Generate matrix of constraint functions
c = zeros(nConstraints*N,1);

% Counter
k = 1; % Column
m = 1; % Row
for i = 1:N
    % System dynamics constraints
    c(m, 1) = x(k+nStates) - x(k) - x(k+10)*x(k+2)*cos(x(k+3));
    c(m+1, 1) = x(k+1+nStates) - x(k+1) - x(k+10)*x(k+2)*sin(x(k+3));
    c(m+2, 1) = x(k+2+nStates) - x(k+2) - x(k+10)*x(k+4);
    c(m+3, 1) = x(k+3+nStates) - x(k+3) - x(k+10)*x(k+5);
    m = m + nConstraintsSystem;
    
    % Constraints stick to triplet (linear interpolation)
        % Safety margin
    if (xInner(i) == xOuter(i) || i == 1)
        c(m,1) = 0;
    else
        c(m,1) = x(k+1) - x(k)*(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i)) - (yInner(i)-(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i))*xInner(i)); 
    end
    m = m + nConstraintsTriplet;
    
    % Obstacle Constraints
    %   (A*p_i - b)^T \lambda > d_{min}
    c(m,1) = x(k+6)*(A(1,1)*x(k) + A(1,2)*x(k+1) - b(1,1))+ x(k+7)*(A(2,1)*x(k) + A(2,2)*x(k+1) - b(2,1)) + x(k+8)*(A(3,1)*x(k) + A(3,2)*x(k+1) - b(3,1)) + x(k+9)*(A(4,1)*x(k) + A(4,2)*x(k+1) - b(4,1));
    %   ||A^T \lambda|| <= 1
    c(m+1,1) = (x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1))^2 + (x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2))^2; 
    m = m + nConstraintsObject;
    
    % Next Interation
    k = k + nStates;
end

% ----------------------------------------------------------------------
function J = jacobian (x,auxdata)

% Define Constants
[N, nStates, ~, ~, ~, ~, xOuter, yOuter, xInner, yInner, A, b] = deal(auxdata{:});

% Number of constraints
nConstraintsSystem = 4;
nConstraintsTriplet = 1;
nConstraintsObject = 2;
nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

% Allocate space for Jacobian Sparse
J = spalloc(nConstraints*N,nStates*(N+1),32*N); 

% Counters
k = 1; % Column
m = 1; % Row

% Fill in Jacobian Matrix
for i = 1:N
    % Jacobian of System Dynamics
    J(m,k) = -1; J(m,k+2) = -x(k+10)*cos(x(k+3)); J(m,k+3) = x(k+10)*x(k+2)*sin(x(k+3)); J(m,k+nStates) = 1; J(m,k+10) = -x(k+2)*cos(x(k+3));
    J(m+1,k+1) = -1; J(m+1,k+2) = -x(k+10)*sin(x(k+3)); J(m+1,k+3) = -x(k+10)*x(k+2)*cos(x(k+3)); J(m+1,k+1+nStates) = 1; J(m+1,k+10) = -x(k+2)*sin(x(k+3)); 
    J(m+2,k+2) = -1; J(m+2,k+4) = -x(k+10); J(m+2,k+2+nStates) = 1; J(m+2,k+10) = -x(k+4);
    J(m+3,k+3) = -1; J(m+3,k+5) = -x(k+10); J(m+3,k+3+nStates) = 1; J(m+3,k+10) = -x(k+5);
    m = m + nConstraintsSystem;

    % Jacobian of Sticking to Triplet Constraint
    if (xInner(i) ~= xOuter(i))
        J(m,k) = -(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i));
        J(m,k+1) = 1;
    end
    m = m + nConstraintsTriplet;
    
    % Jacobian of Obstacle Constraints
    J(m, k) = x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1);
    J(m, k+1) = x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2);
    J(m, k+6) = A(1,1)*x(k) + A(1,2)*x(k+1) - b(1,1);
    J(m, k+7) = A(2,1)*x(k) + A(2,2)*x(k+1) - b(2,1);
    J(m, k+8) = A(3,1)*x(k) + A(3,2)*x(k+1) - b(3,1);
    J(m, k+9) = A(4,1)*x(k) + A(4,2)*x(k+1) - b(4,1);
    
    J(m+1, k+6) = 2*A(1,1)*(x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1)) + 2*A(1,2)*(x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2));
    J(m+1, k+7) = 2*A(2,1)*(x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1)) + 2*A(2,2)*(x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2));
    J(m+1, k+8) = 2*A(3,1)*(x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1)) + 2*A(3,2)*(x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2));
    J(m+1, k+9) = 2*A(4,1)*(x(k+6)*A(1,1) + x(k+7)*A(2,1) + x(k+8)*A(3,1) + x(k+9)*A(4,1)) + 2*A(4,2)*(x(k+6)*A(1,2) + x(k+7)*A(2,2) + x(k+8)*A(3,2) + x(k+9)*A(4,2));
    
    m = m + nConstraintsObject;
    
    % Next Iteration
    k = k + nStates;
end

% ----------------------------------------------------------------------
function J_struct = jacobian_struct(auxdata) 
% Define Constants
[N, nStates] = deal(auxdata{1:2});
nConstraintsTriplet = 1;
nConstraintsSystem = 4;
nConstraintsObject = 2;
nConstraints = nConstraintsTriplet + nConstraintsSystem + nConstraintsObject;

% Generate structure of Jacobian
J = zeros(nConstraints*N,nStates*(N+1));

% Counter
k = 1; % Row
m = 1; % Column

% Fill in Jacobian structure Matrix
for i = 1:N
    % Jacobian of System Dynamics
    J(m,k) = 1; J(m,k+2) = 1; J(m,k+3) = 1; J(m,k+nStates) = 1; J(m,k+10) = 1;
    J(m+1,k+1) = 1; J(m+1,k+2) = 1; J(m+1,k+3) = 1; J(m+1,k+1+nStates) = 1; J(m+1,k+10) = 1;
    J(m+2,k+2) = 1; J(m+2,k+4) = 1; J(m+2,k+2+nStates) = 1; J(m+2,k+10) = 1;
    J(m+3,k+3) = 1; J(m+3,k+5) = 1; J(m+3,k+3+nStates) = 1; J(m+3,k+10) = 1;
    m = m + nConstraintsSystem;   
    
    % Jacobian of Sticking to Triplet Constraint
    J(m,k) = 1;    
    J(m,k+1) = 1;
    m = m + nConstraintsTriplet;
    
    % Jacobian of Obstacle Constraints
    J(m, k) = 1;
    J(m, k+1) = 1;
    J(m, k+6) = 1;
    J(m, k+7) = 1;
    J(m, k+8) = 1;
    J(m, k+9) = 1;
    J(m+1, k+6) = 1;
    J(m+1, k+7) = 1;
    J(m+1, k+8) = 1;
    J(m+1, k+9) = 1; 
    m = m + nConstraintsObject;
    
    % Next Iteration
    k = k + nStates;
end
    
% Convert Jacobian into sparse
    J_struct = sparse(J);