function [xiPred] = saturation(xiPred, N, constr)
%%%% This function is to update xiPred within CONSTRAINER
%%%% Input: Q, R, N, M, T, xiRef, centerPoints, xi, xiPred, constr
%%%%    xiPred: Predicted value of input.
%%%%    N: Length of horizon
%%%%    constr: Constraint construct.
%%%% Output: xiNext
%%%%    xiPred: Next guessed value of input.

%   Ensure that all inputs along horizon are within constraints
    OMEGA_I = 1;
    A_I = 2;

    for i = 1:N % ignore first input since that corresponds to input sent to plant at previous time
        % Check omega
        if xiPred(OMEGA_I,i) > constr.omegaMax
            xiPred(OMEGA_I,i) = constr.omegaMax;
        elseif xiPred(OMEGA_I,i) < constr.omegaMin
            xiPred(OMEGA_I,i) = constr.omegaMin;
        end
        % Check acceleration
        if xiPred(A_I,i) > constr.aMax
            xiPred(A_I,i) = constr.aMax;
        elseif xiPred(A_I,i) < constr.aMin
            xiPred(A_I,i) = constr.aMin;
        end
    end
    
end

