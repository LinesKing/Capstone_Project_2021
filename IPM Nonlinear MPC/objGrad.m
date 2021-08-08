%% obtain objective function gradient
function [g] = objGrad(x,s)

g(1) = x(2);
g(2) = (5+x(1));

ns = size(s,2);
g = [g zeros(1,ns)];