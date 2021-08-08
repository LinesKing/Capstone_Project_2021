function [th] = theta(x,s,bL,bU)
%%%% This function is sum of absolute value for residuals
%%%% Input: x, s, bL, bU
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: th
%%%%    th: time derivative of states

    th = sum(abs(res(x,s,bL,bU)));
end
