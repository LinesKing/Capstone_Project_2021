%% Reference Tracking Race Track
% Using LTV MPC to follow center line path of race track with boundary
% constraints. Simulation exits after all times are simulated.

% Track ('track.mat') sourced from https://github.com/alexliniger/MPCC.
% Default qp solver for MPC is OSQP from https://www.gurobi.com/.

% Unicycle model state variables
% States
% x - position along x-axis
% y - position along y-axis
% theta - heading of unicycle
% v - velocity

% Inputs
% omega - angular velocity of unicycle
% a - acceleration of unicycle

clc;
clear all;
close all;

%% Basic setting
% Load path, outer boundary and center of track.
load('track.mat');

% Track information
xPath = track.center(1,:);
yPath = track.center(2,:);
xCenter = track.center(1,:);
yCenter = track.center(2,:);
xInner = track.inner(1,:);
yInner = track.inner(2,:);
xOuter = track.outer(1,:);
yOuter = track.outer(2,:);
centerPoints = [xCenter; yCenter];
innerPoints = [xInner; yInner];
outerPoints = [xOuter; yOuter];

% Reference
xRef = xPath;
yRef = yPath;
thetaRef = zeros(1, length(xPath));
vRef = zeros(1, length(xPath));
xiRef = [xPath; yPath; thetaRef; vRef];

% Bound constraint of track
xTrackMin = min([xInner, xOuter]);
xTrackMax = max([xInner, xOuter]);
yTrackMin = min([yInner, yOuter]);
yTrackMax = max([yInner, yOuter]);

% Indicate whether to plot vehicle at each time step in the loop. 
% Slows down simulation, but allows debugging of weird behaviour.
flagPlot = true;
 
%% Constraints
% Specify system constraints
constr = struct;

% State constraints
% Note that x-y position constraints cannot be set here since they are 
% constrained by the track boundaries in the controller
constr.thetaMax = inf;  % max heading angle
constr.thetaMin = -inf;  % min heading angle
constr.vMax = 2; % max velocity
constr.vMin = 0; % min velocity

% Input constraints
constr.aMax = 1; % max acceleration
constr.aMin = -constr.aMax; % Min acceleration
constr.omegaMax = 5; % max angular velocity
constr.omegaMin = -constr.omegaMax; % Min angular velocity

%% MPC setting
N = 20;  % length of horizon
M = 1;  % number of iterations to repeat LTV for single time
T = 0.05;  % sampling period
tspan = [0 0];  % loop through times

% Set simulation time
startTime = 0;
endTime = T*(length(xRef)-1);
times = startTime: T: endTime;

% Constant terms
nState = 4;  % number of states in unicycle model
nInput = 2;  % number of inputs in unicycle model
nPosition = 2;  % number of position states in unicycle model
nEta = nState-nPosition; % number of internal states in unicycle model

% Initial states/inputs
vInitial = 0;
omegaInitial = 0;
aInitial = 0;
thetaInitial = atan((yRef(2) - yRef(1))/(xRef(2) - xRef(1)));
xi = [xRef(1), yRef(1), thetaInitial, vInitial]';
u = [omegaInitial, aInitial]';
xiPred = zeros(length(xi), N+1); % guess for states across horizon
uPred = repmat(u, 1, N+1); % guess for inputs across horizon

% Specify Weights
Q = 20*diag([1 1 0 0]); % [x weight, y weight, theta weight, v weight]
R = 20*diag([1 1]); % [omega weight, a weight]

% Test LTV MPC for different horizons
% for n = 1:N
%     [uPred(:,1:n+1), xiPred] = LTVMPC(Q, R, n, M, T, xiRef, ...
%                             centerPoints, xi, uPred(:,1:n+1), constr);
%     % Plot prediction along horizon
%     scatter(xiPred(1,1:n+1), xiPred(2,1:n+1),'+')
%     hold on
%     axis([xTrackMin xTrackMax yTrackMin yTrackMax])
%     plot(xOuter,yOuter,'k')
%     plot(xInner,yInner,'k')
%     hold off
%     pause(1e-1)
% end

%% MPC Simulation
uHis = []; % History of implemented inputs
xiHis = xi; % History of implemented states
epsSim = []; % Record array of epsilons
tspan = [0 0];

% Fill references and center points
xiRef = [xiRef, xiRef(:, 1:N)];
centerPoints = [centerPoints, centerPoints(:, 1:N)];

% fprintf('Simulation Progress:         ')
for k = 1: length(times) % starting off at first time already (not zeroth)
    % fprintf('\b\b\b\b\b\b\b\b%6.2f %%',(i/(length(TIME)-1)*100));
    
    % Convert NMPC to LTV MPC and predict the states and inputs within
    % the horizon
    [uPred, xiPred, epsCurr] = LTVMPC(Q, R, N, M, T, ...
       xiRef(:,k:k+N), centerPoints(:,k:k+N), xi, uPred(:,1:N+1), constr);
    
    % Calculate current states with a ODE model
    tspan = [tspan(2) tspan(2)+T];  % time interval
    [t,xiODE] = ode45(@(t, xiODE) unicycleODE(t,xiODE, ...
                uPred(1:nInput,2)), tspan, xi(1:nState,1));
        
    xi = xiODE(length(xiODE),:)';  % update current states
    xiHis = [xiHis, xi];  % update states history
    uHis = [uHis, uPred(:,2)];  % update inputs history
    epsSim = [epsSim, epsCurr];  % update epsilons
    
    % Plot position states inside simulation
    if flagPlot
        % Plot predicted positions from MPC along horizon
        plot(xiPred(1,1:size(xiPred,2)), ...
                xiPred(2,1:size(xiPred,2)), 'LineWidth', 1.5)
         
        hold on 
        
        % Plot references
        plot(xiRef(1,k+1:k+N),xiRef(2,k+1:k+N))
        
        % Set axis to track edges
        axis([xTrackMin xTrackMax yTrackMin yTrackMax])

        % Plot track boundaries
        plot(xOuter,yOuter,'k')
        plot(xInner,yInner,'k')
        
        % Plot history of positions from t=0 to current t        
        plotTrVec = [xi(1:2); 0];
        plotRot = axang2quat([0 0 1 xi(3)]);
        plotTransforms(plotTrVec', plotRot, "MeshFilePath", ...
                       "groundvehicle.stl", "Parent", gca, "View",...
                       "2D", "FrameSize", 0.15);
        light;
        
        % scatter(xiHis(1,:), xiHis(2,:))
        
        hold off
        
        legend('Predicted pose', 'Reference track')
        
        pause(1e-6)
    end
    
    % Prepare guess of inputs and states along horizon for next iteration
    uPred = [uPred(:,2:N+1), uPred(:,N+1)];
    xiPred = [xi, zeros(length(xi), N)];
end

% Average of final input differences at each time
eps_avg = mean(epsSim);
MSE = immse(xiHis(1:2, 1:length(xRef)), xiRef(1:2, 1:length(xRef)));
fprintf('\n')

%% Plot Positions, States and Inputs
% Plot history of positions during simulation as well as track
figure
plot(xiHis(1,:),xiHis(2,:))
hold on
% axis([track_xmin track_xmax track_ymin track_ymax])
plot(xCenter,yCenter,'r')
plot(xOuter,yOuter,'k')
plot(xInner,yInner,'k')
title('Positions')
xlabel('X')
ylabel('Y')
legend('vehicle','reference', 'boundary')

% Plot history of states during simulation
% Plot x history
figure
subplot(4,1,1)
plot(times(1:end), xiHis(1,1:end-1))
title('States')
ylabel("x [m]")
grid on
% Plot y history
subplot(4,1,2)
plot(times(1:length(times)), xiHis(2,1:end-1))
ylabel("y [m]")
grid on
% Plot v history
subplot(4,1,3)
plot(times(1:length(times)), xiHis(3,1:end-1))
ylabel("\theta (\circ)")
grid on
% Plot theta history
subplot(4,1,4)
plot(times(1:length(times)), xiHis(4,1:end-1))
grid on
ylabel("v [m/s]")

% Plot history of inputs during simulation
figure
% Plot omega history
subplot(2,1,1)
stairs(times(1:length(times)), uHis(1,:))
title('Inputs')
ylabel("\omega [rad/s]")
grid on
% Plot a history
subplot(2,1,2)
stairs(times(1:length(times)), uHis(2,:))
grid on
ylabel("a [m/s^2]")
xlabel("Time (s)")

%% Save data for analyis
record = struct;
record.xiHis = xiHis;
record.uHis = uHis;
record.times = times;
record.constr = constr;
record.track = track;

