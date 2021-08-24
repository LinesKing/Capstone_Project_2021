%% Reference Generating Race Track
% Using OBCA to generate path of partial race track with collision
% avoidance. Simulation exits after all times are simulated.

clc;
clear all;
close all;

%% Basic setting
% Load track a
% load('../tracks/track_a.mat');
% centerPoints = round(track_a.center,4);
% innerPoints = round(track_a.inner,4);
% outerPoints = round(track_a.outer,4);

% Load track b
load('../tracks/track_b.mat');

% Get part of track 
%   (n = 1, m = 40: horizontal line);
%   (n = 40, m = 81: curve);
%   (n = 205, m = 235: vertical line);
n = 1;
m = 40;
centerPoints = round(track_b.center(:,n:m),4);
innerPoints = round(track_b.inner(:,n:m),4);
outerPoints = round(track_b.outer(:,n:m),4);

% Choose starting N point
starting_N_point = 1;

% Receeding Horizon
N = length(centerPoints)-1;
ending_N_point = starting_N_point + N;

% Enable warm start with center points
warm = struct();
warm.xWarm = centerPoints(1,starting_N_point:ending_N_point);
warm.yWarm = centerPoints(2,starting_N_point:ending_N_point);

% Extract inner and outer bounds
% Note: Even though Receeding Horizon is N, we need to feed in N+1 terms, 
% due to initialising all N+1 terms.
xOuter = outerPoints(1,starting_N_point:ending_N_point);
yOuter = outerPoints(2,starting_N_point:ending_N_point);
xInner = innerPoints(1,starting_N_point:ending_N_point);
yInner = innerPoints(2,starting_N_point:ending_N_point);

% Reformulate bounds within margain
safeMargin = 0.05;
for i = 1:N+1
    dist =  norm([xInner(i) yInner(i)]-[xOuter(i) yOuter(i)]);
    t = safeMargin/dist;
    xInnerMargin(i) = (1-t)*xInner(i) + t*xOuter(i);
    yInnerMargin(i) = (1-t)*yInner(i) + t*yOuter(i);
    xOuterMargin(i) = (1-t)*xOuter(i) + t*xInner(i);
    yOuterMargin(i) = (1-t)*yOuter(i) + t*yInner(i);
end

%% Plotting obstacle and track
% Plot track bounds with margin
figure(1)
hold on
margin = plot(xInnerMargin,yInnerMargin, '-.r');
plot(xOuterMargin,yOuterMargin, '-.r');
bounds = plot(xOuter,yOuter, 'k');
plot(xInner,yInner, 'k');

axis([min(min(xInner), min(xOuter)) ...
    max(max(xInner), max(xOuter)) ...
    min(min(yInner), min(yOuter)) ...
    max(max(yInner), max(yOuter))]);
hold on

% axis([-1.5 2 -2 2]) 
% pause()

%% Obstacle constraints
% Obstacle Parameters
% Obstacle 1
origin = [0.4, 0.6]';
theta = deg2rad(0);
length = 0.1/2;
width = 0.15/2;

object = struct();
[A1, b1] = obstacleMatrices(origin, theta, length, width);

% Obstacle 2
origin = [0.9, 0.7]';
theta = deg2rad(0);
length = 0.1/2;
width = 0.15/2;

object = struct();
[A2, b2] = obstacleMatrices(origin, theta, length, width);

% Obstacle 3
origin = [1.3, 0.65]';
theta = deg2rad(0);
length = 0.1/2;
width = 0.15/2;

object = struct();
[A3, b3] = obstacleMatrices(origin, theta, length, width);

% Calculate Ax <= b
object.A = [A1; A2; A3];
object.b = [b1; b2; b3];

% Choose d_min
dMin = 0.05;

%% Plotting obstacle and track
figure(1); 
hold on; grid on;
plotConvexRegion(-A1,-b1,[-1.5 -2],[2 2], 'g', 0.1)
plotConvexRegion(-A2,-b2,[-1.5 -2],[2 2], 'g', 0.1)
plotConvexRegion(-A3,-b3,[-1.5 -2],[2 2], 'g', 0.1)

% Plot track bounds with margin
margin = plot(xInnerMargin,yInnerMargin, '-.r');
plot(xOuterMargin,yOuterMargin, '-.r');
bounds = plot(xOuter,yOuter, 'k');
plot(xInner,yInner, 'k');

axis([min(min(xInner), min(xOuter)) ...
    max(max(xInner), max(xOuter)) ...
    min(min(yInner), min(yOuter)) ...
    max(max(yInner), max(yOuter))]);
hold on

axis([0 1.8 0 1.2]) 
% pause()

%% Solving optimal problem
% Initial Starting Point
initAngle = atan2(warm.yWarm(2)-warm.yWarm(1), warm.xWarm(2)-warm.xWarm(1));
initPoint = [warm.xWarm(1), warm.yWarm(1), initAngle]; 

% Solve optimal path
profile on
solution = OBCAMultiObsTimeOptimalPathPlanning(xInnerMargin, yInnerMargin, xOuterMargin, yOuterMargin, warm, initPoint, object, dMin);
profile viewer
profile off

% Extract Optimal Path
xPath = solution(:, 1);
yPath = solution(:, 2);

%% Plotting optimal track
% Plot path 
figure(1)

optimalPath = scatter(xPath, yPath, 'b*');
hold off
legend([bounds margin optimalPath], 'Track Bounds', 'Safety Margin', 'Optimal Path')
xlabel('x[m]')
ylabel('y[m]')

% plot states and inputs
vPath = solution(:, 3);
thetaPath = solution(:, 4);

aPath = solution(:, 5);
omegaPath = solution(:, 6);

% % Plot history of states during simulation
% figure(2)
% % Plot x path
% subplot(4,1,1)
% plot(xPath)
% title('States')
% ylabel("x [m]")
% grid on
% % Plot y path
% subplot(4,1,2)
% plot(yPath)
% ylabel("y [m]")
% grid on
% % Plot theta
% subplot(4,1,3)
% plot(thetaPath)
% ylabel("\theta (\circ)")
% grid on
% % Plot v
% subplot(4,1,4)
% plot(vPath)
% ylabel("v [m/s]")
% grid on
% 
% % Plot history of inputs during simulation
% figure(3)
% % Plot omega history
% subplot(2,1,1)
% plot(omegaPath)
% ylabel("a [m/s^2]")
% grid on
% % Plot a history
% subplot(2,1,2)
% plot(aPath)
% title('Inputs')
% ylabel("\omega [rad/s]")
% grid on
