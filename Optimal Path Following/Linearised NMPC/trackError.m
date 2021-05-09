function [e, p, j_p] = trackError(N, xiRef, centerPoints, xiPred)
%%%% This function is to calculate error.
%%%% Input: Q, R, N, M, T, xiRef, centerPoints, xi, uGuess, constr
%%%%    N: Length of horizon
%%%%    xiRef: Reference value of state.
%%%%    xiPred: Predicted value of state.
%%%% Output: e, p, j_p
%%%%    e: Error vectors.
%%%%    p: Error vectors in Frenet coordinate.
%%%%    j_p: Transition matrix between Global and Frenet coordinate.

    e = cell(1, size(xiPred, 2));
    p = cell(1, size(xiPred, 2));
    j_p = cell(1, size(xiPred, 2));

    for i = 1:N
        e{i} = xiPred(1:2, i)-xiRef(1:2, i);
        xDiff = xiPred(1:2, i)-centerPoints(:, i);
        angle = atan((centerPoints(2, i+1)-centerPoints(2, i))/...
                            (centerPoints(1, i+1)-centerPoints(1, i)));
        if centerPoints(1, i+1) < centerPoints(1, i)
            angle = angle-pi/2;
        else
            angle = angle+pi/2;
        end
        j_p{i} = [cos(angle) -sin(angle); sin(angle) cos(angle)];
        p{i} = j_p{i}*xDiff;
    end
    
end
