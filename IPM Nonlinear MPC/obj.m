function [ob] = obj(x)
%%%% This function is to achieve the objective function value
%%%% Input: x
%%%%    x - primal variable
%%%% Output: ob
%%%%    ob: object

    ob = x(2)*(5+x(1));
end