function [g] = objGrad(x,s)
%%%% This function is to obtain objective function gradient
%%%% Input: x, s
%%%%    x - primal variable
%%%%    s - slack variable
%%%% Output: r
%%%%    g: gradient of object function

    g(1) = x(2);
    g(2) = (5+x(1));

    ns = size(s,2);
    g = [g zeros(1,ns)];
end