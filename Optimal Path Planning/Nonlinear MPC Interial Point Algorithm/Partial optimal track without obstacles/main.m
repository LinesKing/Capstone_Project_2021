%% Reference Generating Race Track
% Using OBCA to generate optimal path of race track. 
% Simulation exits after all times are simulated.

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
%   (n = 210, m = 235: vertical line);
n = 1;
m = 80;

centerPoints = round(track_b.center(:,n:m),4);
innerPoints = round(track_b.inner(:,n:m),4);
outerPoints = round(track_b.outer(:,n:m),4);

% Choose starting N point
starting_N_point = 1;

% Receeding Horizon
N = length(centerPoints)-1;
ending_N_point = starting_N_point + N;

% Enable warm start with center points
% warm = struct();
% warm.tWarm = 0.002*ones(1,N+1);
% warm.xWarm = centerPoints(1,starting_N_point:ending_N_point) + normrnd(0,0,[1,N+1]);
% warm.yWarm = centerPoints(2,starting_N_point:ending_N_point) + normrnd(0,0.00,[1,N+1]);
% warm.vxWarm = [diff(warm.xWarm) 0]./warm.tWarm;
% warm.vyWarm = [diff(warm.yWarm) 0]./warm.tWarm;
% warm.axWarm = [diff(warm.vxWarm) 0]./warm.tWarm;
% warm.ayWarm = [diff(warm.vyWarm) 0]./warm.tWarm;

% IPOPT solution as warm start
load('../warm_start/IPOPTsol_1to80.mat');
warm.tWarm = solution(starting_N_point:ending_N_point,7)';
warm.xWarm = solution(starting_N_point:ending_N_point,1)' + normrnd(0,0.001,[1,N+1]);
warm.yWarm = solution(starting_N_point:ending_N_point,2)' + normrnd(0,0.001,[1,N+1]);
warm.vxWarm = solution(starting_N_point:ending_N_point,3)';
warm.vyWarm = solution(starting_N_point:ending_N_point,4)';
warm.axWarm = solution(starting_N_point:ending_N_point,5)';
warm.ayWarm = solution(starting_N_point:ending_N_point,6)';

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

% axis([-1.5 2 -2 2]) 
% pause()

%% Solving optimal problem
% Initial Starting Point
initAngle = atan2(warm.yWarm(2)-warm.yWarm(1), warm.xWarm(2)-warm.xWarm(1));
initPoint = [warm.xWarm(1), warm.yWarm(1), initAngle]; 

% Solve optimal path
profile on
[solution, ob_his]= TimeOptimalPathPlanning(xInnerMargin, yInnerMargin, xOuterMargin, yOuterMargin, warm, initPoint);
profile viewer
profile off

% Extract Optimal Path
xPath = solution(:, 1);
yPath = solution(:, 2);

%% Plotting optimal track
% Plot path 
figure(1)

optimalPath = scatter(xPath, yPath, 'b*');
warmPath = scatter(warm.xWarm, warm.yWarm, 'r*');
hold off
% legend([bounds margin optimalPath], 'Track Bounds', 'Safety Margin', 'Optimal Path')
legend([bounds margin optimalPath warmPath], 'Track Bounds', 'Safety Margin', 'Optimal Path', 'Warm Start Path')
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

% Plot objective history
figure(2)
plot(ob_his);

