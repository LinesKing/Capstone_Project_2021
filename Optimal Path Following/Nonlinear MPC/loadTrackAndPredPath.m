%%%% Load tracks and predicted path
load('track.mat');

% Reference
xPath = track.center(1,:);
yPath = track.center(2,:);

% Track bounds
xCenter = track.center(1,:);
yCenter = track.center(2,:);
xInner = track.inner(1,:);
yInner = track.inner(2,:);
xOuter = track.outer(1,:);
yOuter = track.outer(2,:);

load('xpred.mat');  % named xPred
load('ypred.mat');  % named yPred