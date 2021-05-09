function [xiRef] = generateTrackRef(waypoints, initState, termState, ts)
%%%% This function is to fits a smooth, piecewise, continuous curve to a set of waypoints given as [x y] or [x y theta]. 
%%%% After fitting, points along the curve, the path points are expressed as [x y theta kappa dkappa s], where:
%%%%    x y and theta— SE(2) state expressed in global coordinates, with x and y in meters and theta in radians
%%%%    kappa — Curvature, or inverse of the radius, in meters
%%%%    dkappa — Derivative of curvature with respect to arc length in meters per second
%%%%    s — Arc length, or distance along path from path origin, in meters
%%%% Input: waypoints, initState, termState, ts
%%%%    waypoints: Path object, specified as a navPath object.
%%%%    initStates: Initial states as [S ds ddS L dL ddL].
%%%%    termState: terminal states as [S ds ddS L dL ddL].
%%%%    ts: time in seconds
%%%% Output: xiRef
%%%%    xiRef: Interpolated path.

    refPath = referencePathFrenet(waypoints);  % generate a reference path from a set of waypoints
    connector = trajectoryGeneratorFrenet(refPath);  % create a trajectoryGeneratorFrenet object from the reference path.
    [~,trajGlobal] = connect(connector, initState, termState, ts);  % generate a ts second trajectory.
    xiRef = trajGlobal.Trajectory(:,1:3)';
    
end

