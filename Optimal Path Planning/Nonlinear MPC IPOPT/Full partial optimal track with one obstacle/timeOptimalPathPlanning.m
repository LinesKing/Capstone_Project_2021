function [x] = timeOptimalPathPlanning(xInner, yInner, xOuter, yOuter, warm, initPoint, finalPoint)
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

% Define state & input contraints
vxMin = -2;
vyMin = -2;
axMin = -1.5;
ayMin = -1.5;
vxMax = 20;
vyMax = 20;
axMax = 1.5;
ayMax = 1.5;

% Define objetive function weightings
weightDist = 0.1;
weightVel = 0.5;
weightTime = 1;

% Define constants (number of variables)
nSystemStates = 4;
nControlInputs = 2;
nTimeVariable = 1;
nStates = nSystemStates + nControlInputs + nTimeVariable;

% Initialise decision variables, lower bounds and upper bounds
x0 = [];
options.lb = [];
options.ub = [];
for i = 1:N+1
    % Set initial start point to center of track (warm start)
    % [x, y, vx, vy, ax, ay, t]
    x0 = [x0, warm.xWarm(i), warm.yWarm(i), zeros(1, nStates-2)];
   
    % Upper and lower bounds of decision variables (not constraint equations)
    % [x, y, vx, vy, ax, ay]
    if (i == 1) % satisfy initial conditions contraint
        options.lb = [options.lb, [initPoint(1) initPoint(2) vxMin vyMin axMin ayMin]]; %2D array with each row filled with column of vector current states
        options.ub = [options.ub, [initPoint(1) initPoint(2) vxMax vyMax axMax ayMax]];  
    else if (i == N)
            options.lb = [options.lb, [finalPoint(1) finalPoint(2) vxMin vyMin axMin ayMin]]; %2D array with each row filled with column of vector current states
            options.ub = [options.ub, [finalPoint(1) finalPoint(2) vxMax vyMax axMax ayMax]];
        else
            yLb = min(yOuter(i),yInner(i));
            yUb = max(yOuter(i),yInner(i));
            xLb = min(xOuter(i),xInner(i));
            xUb = max(xOuter(i),xInner(i));
            options.lb = [options.lb, [xLb yLb vxMin vyMin axMin ayMin]];
            options.ub = [options.ub, [xUb yUb vxMax vyMax axMax ayMax]]; 
        end
    end
   
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
end

% Set up the auxiliary data
options.auxdata = {N nStates weightDist weightVel weightTime ...
                   xOuter yOuter xInner yInner};

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
[N, nStates, weightDist, weightVel, weightTime] = deal(auxdata{1:5});

% Generate objective function:
%   \Sigma_{i=1}^N (\delta_x_i^T W_d \delta_x_i + 
%   \delta_y_i^T W_d \delta_y_i + \delta_v_x_i^T W_v 
%   \delta_v_x_i + \delta_v_y_i^T W_v \delta_v_y_i + W_t t_i)
f = 0;
k = 1;
for i = 1:N
    f = f + weightDist*(x(k) - x(k+nStates))^2 ...
          + weightDist*(x(k+1) - x(k+1+nStates))^2 ...
          + weightVel*(x(k+2) - x(k+2+nStates))^2 ...
          + weightVel*(x(k+3) - x(k+3+nStates))^2 ...
          + weightTime*x(k+6);
    
    % Next interation
    k = k + nStates;
end

% ----------------------------------------------------------------------
function g = gradient (x,auxdata)
% Define constants
[N, nStates, weightDist, weightVel, weightTime] = deal(auxdata{1:5});

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
end

% ----------------------------------------------------------------------
function c = constraints (x,auxdata)
% Define constants
[N, nStates, ~, ~, ~, xOuter, yOuter, xInner, yInner] = deal(auxdata{:});

nConstraintsSystem = 4;
nConstraintsTriplet = 1;
nConstraints = nConstraintsTriplet + nConstraintsSystem;

% Generate matrix of constraint functions
c = zeros(nConstraints*N,1);

% Counter
k = 1; % Column
m = 1; % Row
for i = 1:N
    % System dynamics constraints
    c(m, 1) = x(k+nStates) - x(k) - x(k+6)*x(k+2);
    c(m+1, 1) = x(k+1+nStates) - x(k+1) - x(k+6)*x(k+3);
    c(m+2, 1) = x(k+2+nStates) - x(k+2) - x(k+6)*x(k+4);
    c(m+3, 1) = x(k+3+nStates) - x(k+3) - x(k+6)*x(k+5);
    m = m + nConstraintsSystem;
    
    % Constraints stick to triplet (linear interpolation)
    if (xInner(i) == xOuter(i) || i == 1)
        c(m,1) = 0;
    else
        c(m,1) = x(k+1) - x(k)*(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i)) - (yInner(i)-(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i))*xInner(i)); 
    end
    m = m + nConstraintsTriplet;
    
    % Next Interation
    k = k + nStates;
end

% ----------------------------------------------------------------------
function J = jacobian (x,auxdata)

% Define Constants
[N, nStates, ~, ~, ~, xOuter, yOuter, xInner, yInner] = deal(auxdata{:});

% Number of constraints
nConstraintsSystem = 4;
nConstraintsTriplet = 1;
nConstraints = nConstraintsTriplet + nConstraintsSystem;

% Allocate space for Jacobian Sparse
J = spalloc(nConstraints*N,nStates*(N+1),N*18); 

% Counters
k = 1; % Column
m = 1; % Row

% Fill in Jacobian Matrix
for i = 1:N
    % Jacobian of System Dynamics
    J(m,k) = -1; J(m,k+2) = -x(k+6); J(m,k+nStates) = 1; J(m,k+6) = -x(k+2);
    J(m+1,k+1) = -1; J(m+1,k+3) = -x(k+6); J(m+1,k+1+nStates) = 1; J(m+1,k+6) = -x(k+3);
    J(m+2,k+2) = -1; J(m+2,k+4) = -x(k+6); J(m+2,k+2+nStates) = 1; J(m+2,k+6) = -x(k+4);
    J(m+3,k+3) = -1; J(m+3,k+5) = -x(k+6); J(m+3,k+3+nStates) = 1; J(m+3,k+6) = -x(k+5);
    m = m + nConstraintsSystem;

    % Jacobian of Sticking to Triplet Constraint
    if (xInner(i) ~= xOuter(i))
        J(m,k) = -(yOuter(i)-yInner(i))/(xOuter(i)-xInner(i));
        J(m,k+1) = 1;
    end
    m = m + nConstraintsTriplet;
    
    % Next Iteration
    k = k + nStates;
end

% ----------------------------------------------------------------------
function J_struct = jacobian_struct(auxdata) 
% Define Constants
[N, nStates] = deal(auxdata{1:2});
nConstraintsTriplet = 1;
nConstraintsSystem = 4;
nConstraints = nConstraintsTriplet + nConstraintsSystem;

% Generate structure of Jacobian
J = zeros(nConstraints*N,nStates*(N+1));

% Counter
k = 1; % Row
m = 1; % Column

% Fill in Jacobian structure Matrix
for i = 1:N
    % Jacobian of System Dynamics
    J(m,k) = 1; J(m,k+2) = 1; J(m,k+nStates) = 1; J(m,k+6) = 1;
    J(m+1,k+1) = 1; J(m+1,k+3) = 1; J(m+1,k+1+nStates) = 1; J(m+1,k+6) = 1;
    J(m+2,k+2) = 1; J(m+2,k+4) = 1; J(m+2,k+2+nStates) = 1; J(m+2,k+6) = 1;
    J(m+3,k+3) = 1; J(m+3,k+5) = 1; J(m+3,k+3+nStates) = 1; J(m+3,k+6) = 1;
    m = m + nConstraintsSystem;   
    
    % Jacobian of Sticking to Triplet Constraint
    J(m,k) = 1;    
    J(m,k+1) = 1;
    m = m + nConstraintsTriplet;
    
    % Next Iteration
    k = k + nStates;
end
    
% Convert Jacobian into sparse
    J_struct = sparse(J);