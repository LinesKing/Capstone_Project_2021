function [controller] = NMPC(N, iterMax, T, Q, R ,constr)
%%%% This function is to generate a nonlinear MPC controller using Yalmip 
%%%% to interface with SNOPT.
%%%% Input: N, iterMax, T, Q, R ,constr
%%%%    N: length of horizon
%%%%    iterMax: maximum number of iterations in IPOPT
%%%%    T: sampling period
%%%%    Q: weight matrix for states
%%%%    R: weight matrix for input (control effort)
%%%%    constr: Constraint construct.
%%%% Output: controller
%%%%    Yalmip object which takes input parameters, solves for decision 
%%%%    variables and returns output parameters

yalmip('clear')

% Parameters
nState = 4; % Number of states
nInput = 2; % Number of inputs
nPosition = 2; % Number of position states

% Decision variables
xi0 = sdpvar(nState,1); % Initial state
u0 = sdpvar(nInput,1); % Initial input
r = sdpvar(repmat(nState,1,N+1),ones(1,N+1)); % Reference signal
e = sdpvar(repmat(nState,1,N+1),ones(1,N+1)); % Reference tracking error
xi = sdpvar(repmat(nState,1,N+1),ones(1,N+1)); % State
u = sdpvar(repmat(nInput,1,N+1),ones(1,N+1)); % Input
%   boundary constraint

% Set up objective function and constraints
constraints = [];

% Constrain initial position
constraints = [constraints, ...
               xi{1} == xi0, ...
               u{1} == u0];

objective = 0;

% Loop over N+1 states
for i = 1:N
    objective = objective + e{i+1}'*Q*e{i+1};
    constraints = [constraints, ... 
        e{i+1} == xi{i+1} - r{i+1}]; % Error between states and reference
   
    % State constraints
    constraints = [constraints, ...
    constr.vMin <= xi{i+1}(4) <= constr.vMax]; % Velocity constraint

end

% Loop over N inputs
for i = 1:N
    objective = objective + (u{i+1}-u{i})'*R*(u{i+1}-u{i});
    % Input constraints
    constraints = [constraints, ...
       constr.omegaMin <= u{i+1}(1) <= constr.omegaMax, ...
       constr.aMin <= u{i+1}(2) <= constr.aMax];
   
    % MPC state constraints
    dx = xi{i}(4)*cos(xi{i+1}(3));
    dy = xi{i}(4)*sin(xi{i+1}(3));
    dtheta = u{i+1}(1);
    dv = u{i+1}(2);
    constraints = [constraints, xi{i+1}(1) == xi{i}(1) + T*dx];
    constraints = [constraints, xi{i+1}(2) == xi{i}(2) + T*dy];
    constraints = [constraints, xi{i+1}(3) == xi{i}(3) + T*dtheta];
    constraints = [constraints, xi{i+1}(4) == xi{i}(4) + T*dv];
end

% Set up controller
parametersIn = {[r{:}], xi0, u0};
solutionsOut = {[u{:}], [xi{:}]};
%ops = sdpsettings('solver','ipopt');
ops = sdpsettings('solver','snopt','usex0',1);
% ops = sdpsettings('solver','fmincon');
controller = optimizer(constraints, objective, ops, ...
                        parametersIn, solutionsOut);
end

