function [Ad, Bd] = unicycleMatrices(xiPred, N, T)
%%%% This function is to generate state and input matrices for each column
%%%% of the array xiPred (each column is new set of states)
%%%% Input:
%%%%    xiPred: 4xN array of states along horizon
%%%%     N: Length of horizon
%%%%     T: Sampling period
%%%% Output:
%%%%    Ad: Discrete state matrices
%%%%    Bd: Discrete input matrices

    Ad = cell(1,N);
    Bd = cell(1,N);
    % Generate state matrices
    for i = 1:N
        Ad{i} = [1, 0, -T*xiPred(4,i)*sin(xiPred(3,i)), T*cos(xiPred(3,i));
                0, 1, T*xiPred(4,i)*cos(xiPred(3,i)), T*sin(xiPred(3,i));
                0, 0, 1, 0;
                0, 0, 0, 1];
    end

    % Generate input matrices
    for i = 1:N
        Bd{i} = [-(T^2*xiPred(4,i)*sin(xiPred(3,i)))/2, (T^2*cos(xiPred(3,i)))/2;
                (T^2*xiPred(4,i)*cos(xiPred(3,i)))/2, (T^2*sin(xiPred(3,i)))/2;
                T, 0;
                0, T];
    end

end

