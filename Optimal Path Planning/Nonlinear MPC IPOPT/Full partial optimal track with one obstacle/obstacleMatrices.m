function [A, b] = obstacleMatrices(origin, theta, length, width)
%%%% This function is to formulate the A and b matrix to represent the 
%%%% obstacle as a compact convex set.
%%%% Input:
%%%%     origin - Center location
%%%%     theta - Rotation degree in rad
%%%%     length - Length
%%%%     width - Width
%%%% Output:
%%%%    A - Matrix
%%%%    b - Vector

A = [-cos(theta) -sin(theta);
      cos(theta)  sin(theta);
      sin(theta) -cos(theta);
     -sin(theta)  cos(theta)];
b = [length/2; length/2; width/2; width/2] + A*origin;
end