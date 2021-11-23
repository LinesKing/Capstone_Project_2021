clc;
clear all;
close all;

% Load track b
load('../tracks/track_b.mat');

% Get part of track 
%   (n = 1, m = 40: horizontal line);
%   (n = 40, m = 81: curve);
%   (n = 210, m = 235: vertical line);
n = 190;
m = 220;

centerPoints = round(track_b.center(:,n:m),4);
innerPoints = round(track_b.inner(:,n:m),4);
outerPoints = round(track_b.outer(:,n:m),4);

% Choose starting N point
starting_N_point = 1;

% Receeding Horizon
% N = length(centerPoints)-1;
N = length(centerPoints)-1;
ending_N_point = starting_N_point + N;

% Enable warm start with center points
warm = struct();
warm.xWarm = centerPoints(1,starting_N_point:ending_N_point);
warm.yWarm = centerPoints(2,starting_N_point:ending_N_point) + normrnd(0,0.001,[1,N+1]);

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

figure(1)
hold on
margin = plot(xInnerMargin,yInnerMargin, '-.r');
plot(xOuterMargin,yOuterMargin, '-.r');
bounds = plot(xOuter,yOuter, '.k');
plot(xInner,yInner, '.k');

axis([min(min(xInner), min(xOuter)) ...
    max(max(xInner), max(xOuter)) ...
    min(min(yInner), min(yOuter)) ...
    max(max(yInner), max(yOuter))]);
hold on

plot(xCenter,yCenter, '.b');

for i = 1:N+1
    plot([xInnerMargin(i),xOuterMargin(i)],[yInnerMargin(i),yOuterMargin(i)],'-g');
end


legend('Safety Margin 1','Safety Margin 2', 'Boundary 1', 'Boundary 2', 'Center point', 'Feasible region')
