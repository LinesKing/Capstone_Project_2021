function [theta] = residualAbsSum(x,s,bL,bU,auxdata)
%%%% This function is sum of absolute value for residuals
%%%% Input: x, s, bL, bU
%%%%    x - primal variable
%%%%    s - slack variable
%%%%    bL: lower bound for vector b 
%%%%    bU: upper bound for vector b
%%%% Output: th
%%%%    th: time derivative of states

    theta = sum(abs(constaintResidual(x,s,bL,bU,auxdata)));
end
