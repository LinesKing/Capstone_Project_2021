function [uPred, xiPred, epsCurr] = LTVMPC(Q, R, N, M, T, xiRef, ...
                                        centerPoints, xi, uPred, constr)
%%%% This function is to control a LTV system by MPC using OSQP.
%%%% Input: Q, R, N, M, T, xiRef, centerPoints, xi, uPred, constr
%%%%    Q: States weight.
%%%%    R: Input weight.
%%%%    N: Length of horizon
%%%%    M: Number of iterations to repeat LTV for single time.
%%%%    xiRef: Reference value of state.
%%%%    centerPoints: Center points.
%%%%    xi: state.
%%%%    uPred: Predicted value of input.
%%%%    constr: Constraint construct.
%%%% Output: xiNext
%%%%    uPred: Next predicteded value of input.
%%%%    xiPred: Next predicteded value of state.
%%%%    epsCurr: Current epsilon over iterations

%   Compute control inputs and states over horizon
    nState = 4;  % number of states in unicycle model
    nInput = 2;  % number of inputs in unicycle model
    nPosition = 2;  % number of position states in unicycle model
    nEta = nState-nPosition; % number of internal states in unicycle model

    omegaGuessCurrStepHis = [uPred(1,1:N+1); zeros(M+1,N+1)];
    aGuessCurrStepHis = [uPred(2,1:N+1); zeros(M+1,N+1)];
    
    for j = 1:M
        if M > 1
            % Take average of inputs ((M * N) as (1 * N)) inside iteration 
            % to be next guess for uPred
            uPred(1, 1:N+1) = mean(omegaGuessCurrStepHis(1:j,:), 1);
            uPred(2, 1:N+1) = mean(aGuessCurrStepHis(1:j,:), 1);
        end
        
        % Calculate predicted states
        xiPred = [xi, zeros(nState, N)];
        for i = 1:N
            xiPred(:,i+1) = LTVPred(xiPred(:,i), uPred(:,i+1), T);
        end
        
        % Get state and input matrices
        [Ad, Bd] = unicycleMatrices(xiPred(:,1:N), N, T);
        
        etaMin = [constr.thetaMin; constr.vMin];
        etaMax = [constr.thetaMax; constr.vMax];
        uMin = [constr.omegaMin; constr.aMin];
        uMax = [constr.omegaMax; constr.aMax];
                
        [e, p, j_p] = trackError(N, xiRef(:,2:N+1), centerPoints(:,1:N+1), xiPred(:,2:N+1));
        
        [P_qp, q_qp, Aineq_qp,l_qp, u_qp, Aeq_qp, c_qp] = ...
            LTVMPC2QP(N, Q, R, nState,nPosition, nEta, nInput,  ...
            uPred(:,1:N+1), xiPred(:,1:N+1), e, p, j_p,...
            etaMin, etaMax, uMin, uMax, Ad, Bd);

        [delta_xiPred, delta_uPred] = QP2OSQP(P_qp, q_qp, Aineq_qp, ...
                l_qp, u_qp, Aeq_qp, c_qp, N, nState, nInput);
        
        % Update with solutions, retain same inputs if solution is unfeasible
        uPred(:,2:N+1) = uPred(:,2:N+1) + delta_uPred(:,2:N+1);
        xiPred(:,2:N+1) = xiPred(:,2:N+1) + delta_xiPred(:,2:N+1);
        [uPred(:,2:N+1)] = saturation(uPred(:,2:N+1), N, constr);
        
        % Record epsilon over iterations
        epsCurr = norm([uPred(1,1:N+1) - omegaGuessCurrStepHis(j,1:N+1);
                         uPred(2,1:N+1) - aGuessCurrStepHis(j,1:N+1)], 2);

        % Record U over iterations
        omegaGuessCurrStepHis(j+1,:) = uPred(1,1:N+1);
        aGuessCurrStepHis(j+1,:) = uPred(2,1:N+1);

    end
end

