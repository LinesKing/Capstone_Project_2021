function [xiPred] = saturation(xiPred, N, constr)
%%%% This function is to update xiPred within CONSTRAINER
%%%% Input: Q, R, N, M, T, xiRef, centerPoints, xi, xiPred, constr
%%%%    xiPred: Predicted value of input.
%%%%    N: Length of horizon
%%%%    constr: Constraint construct.
%%%% Output: xiNext
%%%%    xiPred: Next guessed value of input.

    for i = 1: N % ignore first input since that corresponds to input sent to plant at previous time
        % Check omega
        if xiPred(1, i) > constr.omegaMax
            xiPred(1, i) = constr.omegaMax;
        elseif xiPred(1, i) < constr.omegaMin
            xiPred(1, i) = constr.omegaMin;
        end
        % Check acceleration
        if xiPred(2, i) > constr.aMax
            xiPred(2, i) = constr.aMax;
        elseif xiPred(2, i) < constr.aMin
            xiPred(2, i) = constr.aMin;
        end
    end
    
end

