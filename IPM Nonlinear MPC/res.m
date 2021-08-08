% residual with equality constants and inequality slack variables
function [r] = res(x,s,bL,bU)

% get base residual
r(1) = x(1)*x(2);
r(2) = x(1)^2 + x(2)^2;
j = 0;
for i = 1:size(r,2)
    if (bU(i)==bL(i))
        % equality constant
        r(i) = r(i) - bL(i);
    else
        % inequality slack
        j = j + 1;
        r(i) = r(i) - s(j);
    end
end

