function [Pmat, qmat, Aineq, l, u, Aeq, c] = LTVMPC2QP(N, Q, R, nState, ...
                    nInput, nPosition, nEta, uPred, xiPred, ...
                    e, p, j_p, etaMin, etaMax, uMin, uMax, Ad, Bd)
%%%% This function is to convert LTV MPC into QP problem. 
%%%% Input: N, Q, R, nState, nInput, nPosition, nEta, uPred, xiPred,
%%%%        e, p, j_p, etaMin, etaMax, uMin, uMax, Ad, Bd
%%%%    N: Length of horizon.
%%%%    Q: States weight.
%%%%    R: Input weight.
%%%%    nState: Number of states in unicycle model.
%%%%    nInput: Number of inputs in unicycle model.
%%%%    nPosition: Number of position states in unicycle model.
%%%%    nEta: Number of internal states in unicycle model.
%%%%    uPred: N_INPUTS x N+1 matrix of guesses for inputs from k-1 to 
%%%%            k+N-1.
%%%%    xiPred: N_STATES X N+1 matrix of guesses for inputs from k to k+N.
%%%%    e: N cells of errors with dimension 2 x 1 from k+1 to k+N.
%%%%    p: Error vectors in Frenet coordinate.
%%%%    j_p: Transition matrix between Global and Frenet coordinate.
%%%%    etaMin: min value of internal state.
%%%%    etaMax: max value of internal state.
%%%%    uMin: min value of input.
%%%%    uMax: max value of input.
%%%%    Ad: N cells of state matrices from k to k+N-1
%%%%    Bd: N cells of state matrices from k to k+N-1
%%%% Output: Pmat, qmat, Aineq, l, u, Aeq, c
%%%%    Pmat: Pmat in QP
%%%%    qmat: Pmat in QP
%%%%    Aineq: Aineq in QP
%%%%    l: l in QP
%%%%    u: u in QP
%%%%    Aeq: Aeq in QP
%%%%    c: c in QP

    J_err = [1 0 0 0; 0 1 0 0];
    % Formulate P matrix
    Pmat = zeros((N+1)*nState + (N+1)*nInput, (N+1)*nState + (N+1)*nInput);
    % State section of P matrix
    for i = 1:N
        Pmat(i*nState+1:(i+1)*nState, i*nState+1:(i+1)*nState) = J_err' * Q * J_err;
    end
    % Input section of P matrix
    I1blk = [-eye(nInput); eye(nInput)];
    I1 = kron(eye(N), I1blk);
    I2blk = [eye(nInput); eye(nInput)];
    I2 = blkdiag(eye(nInput), kron(eye(N-1), I2blk), eye(nInput));
    Rblkdiag = kron(eye(N), R);
    Pmat((N+1)*nState+1:end, (N+1)*nState+1:end) = I2'*I1*Rblkdiag*I1'*I2;
    Pmat = 2*Pmat;

    % Formulate q matrix
    qmat = zeros((N+1)*nState + (N+1)*nInput,1);
    % State section of q matrix
    for i = 1:N
        qmat(i*nState+1:(i+1)*nState,1) = 2*J_err'*Q*e{i};
    end
    % Input section of q matrix
    for i = 1:N+1
        if i == 1
            qmat((N+1)*nState+1:(N+1)*nState+nInput,1) = 2*(-R*(uPred(:,2) - uPred(:,1)));
        elseif i == N+1
            qmat((N+1)*nState + (i-1)*nInput + 1:(N+1)*nState + (i)*nInput,1) = 2*(-R*(uPred(:,N) - uPred(:,N+1)));
        else
            qmat((N+1)*nState + (i-1)*nInput + 1:(N+1)*nState + (i)*nInput,1) = 2*(-R*(uPred(:,i+1) - 2*uPred(:,i) + uPred(:,i-1)));
        end
    end


    % Formulate Aeq
    Aeq = zeros(nState + nInput + N*nState, (N+1)*nState + (N+1)*nInput);
    % Equality constraint for initial state
    Aeq(1:nState,1:nState) = eye(nState);
    % Equality constraint for initial input
    Aeq(nState+1:nState+nInput, (N+1)*nState+1:(N+1)*nState+nInput) = eye(nInput);
    % Equality constraints for system dynamics
    for i = 0:N-1
        % State matrix component of constraint
        Aeq(nState+nInput+i*nState+1:nState+nInput+(i+1)*nState, i*nState+1:(i+1)*nState) = Ad{i+1};
        % Input matrix component of constraint
        Aeq(nState+nInput+i*nState+1:nState+nInput+(i+1)*nState, (N+1)*nState+(i+1)*nInput+1:(N+1)*nState+(i+2)*nInput) = Bd{i+1};
        % Identity component of constraint
        Aeq(nState+nInput+i*nState+1:nState+nInput+(i+1)*nState, (i+1)*nState+1:(i+2)*nState) = -eye(nState);
    end

    % Formulate equality bound
    c = zeros(nState + nInput + N*nState, 1);

    % Formulate Aineq
    Aineq = zeros(N*(nState) + N*(nInput), (N+1)*(nState) + (N+1)*(nInput));
    % State section of Aineq matrix

    for i = 1:N
        Aineq((i-1)*nState+1:i*nState, i*nState+1:(i+1)*nState,1) = blkdiag(j_p{i},eye(2));  % [J_p; A_eta; A_beta];
    end
    % Input section of Aineq matrix
    Aineq(N*nState+1:end, (N+1)*nState+nInput+1:end) = kron(eye(N), eye(nInput));

    % Formulate lower inequality bound
    l = zeros(N*nState + N*nInput,1);
    % state component of lower inequality bound
    for i = 1:N
        l((i-1)*nState+1:(i-1)*nState+nPosition,1) = [-0.18 -0.1]'-p{i};
        l((i-1)*nState+nPosition+1:(i-1)*nState+nPosition+nEta,1) = etaMin - xiPred(nPosition+1:nPosition+nEta,i+1);
    end
    % input component of lower inequality bound
    for i = 0:N-1
        l(N*nState+i*nInput+1:N*nState+i*nInput+nInput,1) = uMin - uPred(1:nInput,i+2);
    end

    % Formulate upper inequality bound
    u = zeros(N*nState + N*nInput,1);
    % state component of upper inequality bound
    for i = 1:N
        u((i-1)*nState+1:(i-1)*nState+nPosition,1) = [0.18 0.1]'-p{i} ;
        u((i-1)*nState+nPosition+1:(i-1)*nState+nPosition+nEta,1) = etaMax - xiPred(nPosition+1:nPosition+nEta,i+1);
    end
    % input component of upper inequality bound
    for i = 0:N-1
        u(N*nState+i*nInput+1:N*nState+i*nInput+nInput,1) = uMax - uPred(1:nInput,i+2);
    end
end

