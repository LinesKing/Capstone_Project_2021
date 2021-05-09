function [delta_xi, delta_u] = QP2OSQP(Pmat, qmat, Aineq, l, u, Aeq, c, ...
                                        N, nState, nInput)
%%%% This function is to convert QP matrix to quadprog solvable form and solve.
%%%% Input: Pmat, qmat, Aineq, l, u, Aeq, c, N, nState, nInput
%%%%    Pmat: Pmat in QP
%%%%    qmat: Pmat in QP
%%%%    Aineq: Aineq in QP
%%%%    l: l in QP
%%%%    u: u in QP
%%%%    Aeq: Aeq in QP
%%%%    c: c in QP
%%%%    nState: Number of states in unicycle model.
%%%%    nInput: Number of inputs in unicycle model.
%%%% Output: delta_u, delta_xi
%%%%    delta_u: delta input
%%%%    delta_xi: delta state
    
    % OSQP constraints
    A_OSQP = [Aineq; Aeq];
    l_OSQP = [l; c];
    u_OSQP = [u; c];
    
    % Create an OSQP object
    prob = osqp;
    
    % Setup workspace
    prob.setup(Pmat, qmat, A_OSQP, l_OSQP, u_OSQP, 'warm_start', true, 'verbose', false);
    
    % Solve
    res = prob.solve();
    
     % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end
    
    % Apply first control input to the plant
    delta_xiCol = res.x(1:(N+1)*nState);
    delta_xi = reshape(delta_xiCol, [nState, N+1]);
    delta_uCol = res.x((N+1)*nState+1:(N+1)*nState+(N+1)*nInput);
    delta_u = reshape(delta_uCol, [nInput N+1]);
    
end

