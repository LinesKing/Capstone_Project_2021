%% Reference Generating Race Track
% Using OBCA to generate path of partial race track with collision
% avoidance. Simulation exits after all times are simulated.

clc;
clear all;
close all;

%% Obstacle constraints
% Obstacle Parameters
% Obstacle 1 for horizontal line
% origin = [0.6, 0.64]';
% theta = deg2rad(0);
% length = 0.2;
% width = 0.1;

% Obstacle 2 for curve
origin = [2.3, 0.3]';
theta = deg2rad(60);
length = 0.2;
width = 0.15;

% Obstacle 3 for vertical line
% origin = [-1.35, -0.3]';
% theta = deg2rad(0);
% length = 0.15;
% width = 0.2;

% Calculate Ax <= b
object = struct();
[A, b] = obstacleMatrices(origin, theta, length, width);
object.A = A;
object.b = b;

% Choose d_min
dMin = 0.05;

%% Plotting obstacle
figure(1); 
hold on; grid on;
plotConvexRegion(-A,-b,[-1.5 -2],[4 2], 'r', 0.5)
hold on
plot(origin(1),origin(2),'*')

xlabel('x')
ylabel('y')
