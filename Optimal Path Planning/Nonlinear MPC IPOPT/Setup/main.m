%% Reference Generating Race Track
% Using OBCA to generate optimal path of race track. 
% Simulation exits after all times are simulated.

clc;
clear all;
close all;

%% Basic setting
% Load track a
load('../tracks/track_a.mat');
centerPoints = round(track_a.center,4);
innerPoints = round(track_a.inner,4);
outerPoints = round(track_a.outer,4);

% Load track b
% load('../tracks/track_b.mat');
% centerPoints = round(track_b.center,4);
% innerPoints = round(track_b.inner,4);
% outerPoints = round(track_b.outer,4);

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
xCenter = centerPoints(1,starting_N_point:ending_N_point);
yCenter = centerPoints(2,starting_N_point:ending_N_point);


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
% margin = plot(xInnerMargin,yInnerMargin, '-.r');
% plot(xOuterMargin,yOuterMargin, '-.r');
bounds = plot(xOuter,yOuter, '.k');
plot(xCenter,yCenter, '.b');
plot(xInner,yInner, '.k');

axis([min(min(xInner), min(xOuter)) ...
    max(max(xInner), max(xOuter)) ...
    min(min(yInner), min(yOuter)) ...
    max(max(yInner), max(yOuter))]);
hold on

% axis([-1.5 2 -2 2]) 
% pause()

