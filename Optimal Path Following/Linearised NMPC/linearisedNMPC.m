clc;
clear all;
close all;

%% Basic setting
% Load path, outer boundary and center of track as xPath, yPath, yCenter, 
% xInner, yInner, xOuter, yOuter
loadTrackAndPredPath;  
xRef = xPath;
yRef = yPath;

thetaRef = zeros(1, length(xPath));
centerPoints = [xCenter; yCenter];
innerPoints = [xInner; yInner];
outerPoints = [xOuter; yOuter];
xiRef = [xPath; yPath];
% Edge of track
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
% Specify Weights
Q = 20*diag([1 1]); % [x weight, y weight, theta weight, v weight]
R = 5*diag([1 1]); % [omega weight, a weight]

%% MPC setting
N = 16;  % length of horizon
M = 1;  % number of iterations to repeat LTV for single time
T = 0.05;  % aampling period
tspan = [0 0];  % loop through times

% Set simulation time
startTime = 0;
endTime = T*(length(xRef)-N-1);
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

% Test MPC
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
epsSim = []; % Record array of epsilons, corresponding to norm of 

% Loop through times for simulation
tspan = [0 0];
% fprintf('Simulation Progress:         ')
% difference between inputs in current and previous iteration. Indicates
% the convergence of inputs in an iteration.
for k = 1:length(times) % starting off at first time already (not zeroth)
    % fprintf('\b\b\b\b\b\b\b\b%6.2f %%',(i/(length(TIME)-1)*100));
    [uPred, xiPred, epsCurr] = LTVMPC(Q, R, N, M, T, ...
       xiRef(:,k:k+N), centerPoints(:,k:k+N), xi, uPred(:,1:N+1), constr);

    epsSim = [epsSim epsCurr];

    % Apply input to plant
    tspan = [tspan(2) tspan(2)+T];
    [t,xiODE] = ode45(@(t, xiODE) unicycleODE(t,xiODE, ...
                uPred(1:nInput,2)), tspan, xi(1:nState,1));
    xi = xiODE(length(xiODE),:)';
    
    % Update history
    uHis = [uHis, uPred(:,2)];
    xiHis = [xiHis,xi];
    
    % Plot position inside simulation
    if flagPlot
        % Plot solved positions from MPC along horizon
        plot(xiPred(1,1:size(xiPred,2)), ...
             xiPred(2,1:size(xiPred,2)), 'LineWidth', 1.5)
        hold on 

        % Set axis for plot to track edges
        axis([xTrackMin xTrackMax yTrackMin yTrackMax])

        % Plot track boundaries
        plot(xOuter,yOuter,'k')
        plot(xInner,yInner,'k')
        % Plot references
        plot(xiRef(1,k+1:k+N),xiRef(2,k+1:k+N))

        % Plot history of positions from t=0 to current t
        % scatter(xiHis(1,:),xiHis(2,:))
        
        plotTrVec = [xi(1:2); 0];
        plotRot = axang2quat([0 0 1 xi(3)]);
        plotTransforms(plotTrVec', plotRot, "MeshFilePath", "groundvehicle.stl", "Parent", gca, "View","2D", "FrameSize", 0.15);
        light;
        
        hold off
        % legend('Reference track', 'Current pose')
        pause(1e-6)
    end
    
    % Prepare guess of inputs and states along horizon for next iteration
    uPred = [uPred(:,2:N+1), uPred(:,N+1)];
    xiPred = [xi, zeros(length(xi), N)];
end

% Average of final input differences at each time
eps_avg = mean(epsSim);
fprintf('\n')

%% Plot System, States and Inputs
plotPositions
plotStates
plotInputs

%% Save data for analyis
record = struct;
record.xiHis = xiHis;
record.uHis = uHis;
record.times = times;
record.constr = constr;



