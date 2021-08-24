function [chi] = TimeOptimalPathPlanning(xInner, yInner, xOuter, yOuter, warm, initPoint)
%%%% This function is to find the optimal path by Large-scale Interior Point solution method.
%%%% Input:
%%%% 	xInner - Inner bound of x
%%%%    yInner - Inner bound of y
%%%%    xOuter - Outer bound of x
%%%%    yOuter - Outer bound of y
%%%%    warm - warm starting points
%%%%    initPoint - Initial points
%%%% Output:
%%%%    sol - Optimal solution

    %% Basic settings for problem
    prob = struct;
    prob.nu = 1;  % initial merit function weighting of constraint residuals
    prob.gamma = 0.5;
    prob.mu = 1;  % the initial barrier term
    prob.eta = 0.1;
    prob.epsilon = 1e-8;

    % maximum iterations
    prob.maxIter = 1000;

    prob.tauMax = 0.01; % update tau
    prob.eTol = 1e-3;  % solution error tolerance
    prob.k_mu  = 0.2; % check for new barrier problem (0,1)
    
    % Receeding horizon
    N = length(xOuter) - 1;

    % Max time between points (Must be chosen to something physcally poissble)
    tMax = 0.2; 
    tMin = 1e-4;

    % Define state & input contraints
    vxMin = -2;
    vyMin = -2;
    axMin = -1.5;
    ayMin = -1.5;
    vxMax = 20;
    vyMax = 20;
    axMax = 1.5;
    ayMax = 1.5;

    % Define objetive function weightings
    weightDist = 0.1; % 0.1
    weightVel = 0.5; % 0.5
    weightTime = 1; % 1

    % Define constants (number of variables)
    nSystemStates = 4;
    nControlInputs = 2;
    nTimeVariable = 1;
    nStates = nSystemStates + nControlInputs + nTimeVariable;

    % Set up the auxiliary data
    auxdata = {N nStates weightDist weightVel weightTime xOuter yOuter xInner yInner};
               
    %% Primal dual interior method
    % 1. Initialise decision variables with lower bounds and upper bounds
    chi = [];
    chiL = [];
    chiU = [];
    for i = 1:N+1
        % Set initial start point to center of track (warm start)
        % chi = [x, y, vx, vy, ax, ay, t]
        chi = [chi, warm.xWarm(i), warm.yWarm(i), warm.vxWarm(i), warm.vyWarm(i), warm.axWarm(i), warm.ayWarm(i), warm.tWarm];
        % Upper and lower bounds of decision variables (box constraints)
        if (i == 1) % satisfy initial conditions contraint
            chiL = [chiL, [initPoint(1) initPoint(2) vxMin vyMin axMin ayMin]];
            chiU = [chiU, [initPoint(1) initPoint(2) vxMax vyMax axMax ayMax]];  
        else
            yLb = min(yOuter(i),yInner(i));
            yUb = max(yOuter(i),yInner(i));
            xLb = min(xOuter(i),xInner(i));
            xUb = max(xOuter(i),xInner(i));
            chiL = [chiL, [xLb yLb vxMin vyMin axMin ayMin]];
            chiU = [chiU, [xUb yUb vxMax vyMax axMax ayMax]]; 
        end
        % Lower and Upper Bounds for sampling time
        chiL = [chiL, tMin];
        chiU = [chiU, tMax];
    end
    
    
    % 2. Initialise upper and lower bounds for constraints
    bL = [];
    bU = [];
    for i = 1:N
        % System constraints
        bL = [bL, [0 0 0 0]];
        bU = [bU, [0 0 0 0]];
        % Triplet contraints
        bL = [bL, 0];
        bU = [bU, 0];
    end
    
    
    % 3. Initialise slack variables for inequality constraints
    s = [];  % slack variables
    sL = [];  % lower bound for s
    sU = [];  % upper bound for s
    k = 0;
    % Initialize constaint residuals
    c = constaintResidual(chi,s,bL,bU,auxdata);
    % number of constraints (inequality or equality)
    nc = max(size(c));  
    for i = 1:nc
        if(bU(i)>bL(i)) % no slack variable when bU(i) == bL(i)
            % inequality constraints
            k = k + 1;
            s(k) = c(i);
            sL(k) = bL(i);
            sU(k) = bU(i);
        end
    end

    % Problem size
    nchi = max(size(chi));  % number of primal variables chi
    ns = max(size(s));  % number of slack variables s
   
    % Check and make the problem feasible
    % move chi into feasible region and off boundary initially
    for i = 1:nchi
        if (chi(i)<chiL(i) || chi(i)>chiU(i))
            fprintf(1,'Moving initial chi to interior region\n')
            break
        end
    end
    % move primal variables to be feasible
    for i = 1:nchi
        if (chi(i)<chiL(i))
            chi(i) = (chiU(i)+chiL(i))/2;
        end
        if (chi(i)>chiU(i))
            chi(i) = (chiU(i)+chiL(i))/2;
        end
    end
    % move slack variables to be feasible
    for i = 1:ns
        if (s(i)<sL(i))
            s(i) = (sU(i)+sL(i))/2;
        end
        if (s(i)>sU(i))
            s(i) = (sU(i)+sL(i))/2;
        end
    end
    
    % 4. Initialise dual variables (variable constraint multipliers z)
    % zL = mu / (x-xL)
    for i = 1:nchi
        zL(i) = prob.mu / (chi(i)-chiL(i)+prob.epsilon);
    end
    for i = 1:ns
        zL(nchi+i) = prob.mu / (s(i)-sL(i)+prob.epsilon);
    end
    % zU = mu / (xU-x)
    for i = 1:nchi
        zU(i) = prob.mu / (chiU(i)-chi(i)+prob.epsilon);
    end
    for i = 1:ns
        zU(nchi+i) = prob.mu / (sU(i)-s(i)+prob.epsilon);
    end
    % Initialize objective gradient and constraint jacobian
    g = objectiveGradient(chi,s,auxdata);
    J = constraintJacobian(chi,bL,bU,auxdata);
    % Initialise dual variables (equality constraint multipliers lambda)
    lambda = lsqminnorm(full(J*J'),J)*(zL'-zU'-g');
    lambda = lambda';

    
    % 5. Prepare for iterations
    % Initial alpha for line search
    alphaPr = 0.06;  % primal alpha
    alphaDu = 0.06;  % dual alpha
    alphaPrMax = 0.06;
    alphaDuMax = 0.06;
    % Initialize iteration count
    iter = 0;
    % Print iteration zero
    iprint(prob,iter,chi,lambda,zL,zU,alphaPr,alphaDu,s,chiL,chiU,bL,bU,auxdata);

    
    % 6. Start iterating
    for iter = 1:prob.maxIter

        % Objective function with barrier terms
        phi = barrierObjective(prob,chi,chiL,chiU,s,bL,bU,auxdata);

        % Make diagonal matrices 
        % Dual variables (variable constraint multipliers z)
        ZL = diag(zL);
        ZU = diag(zU);
        % Sigmas
        dL = [chi-chiL+10-8 s-sL+10-8];
        dU = [chiU-chi+10-8 sU-s+10-8];
        DL = diag(dL);
        DU = diag(dU);
        invDL = zeros(nchi+ns,nchi+ns);
        invDU = zeros(nchi+ns,nchi+ns);
        SigmaL = zeros(nchi+ns,nchi+ns);
        SigmaU = zeros(nchi+ns,nchi+ns);
        for i = 1:nchi+ns
           invDL(i,i) = 1 / DL(i,i); 
           invDU(i,i) = 1 / DU(i,i); 
           SigmaL(i,i) = ZL(i,i) / DL(i,i);
           SigmaU(i,i) = ZU(i,i) / DU(i,i);
        end

        % Construct and solve Ax=b
        % Calculate 1st and 2nd derivatives
        g = objectiveGradient(chi,s,auxdata);
        J = constraintJacobian(chi,bL,bU,auxdata);
        W = lagrangianHessian(chi,s,lambda,auxdata);
        H = W + SigmaL + SigmaU;       
        c = constaintResidual(chi,s,bL,bU,auxdata);
        % ones vector
        e = ones(nchi+ns,1);
        % construct A and b
        A = [H,J';J,zeros(nc,nc)];
        b(1:nchi+ns,1) = g' + J'*lambda' - prob.mu*invDL*e + prob.mu*invDU*e;
        b(nchi+ns+1:nchi+ns+nc,1) = c(1:nc)';
        % compute search direction, solving Ax = b
        d = lsqminnorm(A,-b);
        dchi = d(1:nchi);
        for i = 1:nchi
            if chiL(i) == chiU(i)
                dchi(i) = 0;
            end
        end
        ds = d(nchi+1:nchi+ns,1);
        dlambda = d(nchi+ns+1:nchi+ns+nc,1);
        dzL = prob.mu*invDL*e - zL' - SigmaL*[dchi; ds];
        dzU = prob.mu*invDU*e - zU' + SigmaU*[dchi; ds];
        search = [d;dzL;dzU];
    
        % Compute new points
        chiNew = chi + alphaPr * dchi';
        if (ns>=1)
            sNew = s + alphaPr * ds';
        else
            sNew = [];
        end
        zLNew = zL + alphaDu * dzL';
        zUNew = zU + alphaDu * dzU';
        
        % Check constraint violations for alpha
        % Max alpha is that which brings the search point to within "tau" of constraint
        % tau is 0 to 0.01 (tau = mu when mu<0.01, otherwise tau=0.01)
        tau = min(prob.tauMax,100*prob.mu);

        for i = 1:nchi
            if(chiNew(i)<chiL(i))
                alphaPrMax = min(alphaPrMax,(chiL(i)+tau*(chi(i)-chiL(i))-chi(i))/dchi(i,1));
            end
            if(chiNew(i)>chiU(i))
                alphaPrMax = min(alphaPrMax,(chiU(i)-tau*(chiU(i)-chi(i))-chi(i))/dchi(i,1));
            end
        end
        for i = 1:ns
            if(sNew(i)<sL(i))
                alphaPrMax = min(alphaPrMax,(sL(i)+tau*(s(i)-sL(i))-s(i))/ds(i,1));
            end
            if(sNew(i)>sU(i))
                alphaPrMax = min(alphaPrMax,(sU(i)+tau*(sU(i)-s(i))-s(i))/ds(i,1));
            end
        end
        for i = 1:nchi       
            if(zLNew(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zL(i)-zL(i))/dzL(i,1));
                zLNew(i) = 0;
            end
            if(zUNew(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zU(i)-zU(i))/dzU(i,1));
                zUNew(i) = 0;
            end
        end


        % Line search
        % Set alpha as approach to constraint
        alphaDu = min(alphaDu,alphaDuMax);
        alphaPr = min(alphaPr,alphaPrMax);
            
        % Predicted and actual reduction in the merit function
        pred = -alphaPr*g*[dchi;ds] - prob.gamma*alphaPr^2*[dchi;ds]'*H*[dchi;ds] + prob.nu*(norm(c',1)-norm(c'+alphaPr*J*[dchi;ds],1));
        ared = merit(prob,chi,chiL,chiU,s,bL,bU,auxdata) - merit(prob,chiNew,chiL,chiU,sNew,bL,bU,auxdata);  
        
        % Compare actual reduction to predicted reduction as long as the 
        % actual reduction is a fraction of the predicted reduction then
        % accept the trial point.

        while ared < prob.eta*pred 
            % reject point and move alpha_x
            alphaPr = alphaPr * 0.5;
            alphaDu = alphaDu * 0.5;
            % compute new points
            chiNew = chi + alphaPr * dchi';
            if (ns>=1)
                sNew = s + alphaPr * ds';
            else
                sNew = [];
            end
            
            % determine nu
            predPhi1 = g*[dchi;ds];
            % set 2nd derivative contribution to zero if < 0
            predPhi2 = max(0,0.5*[dchi;ds]'*H*[dchi;ds]);
            predPhiDecrease = predPhi1 + predPhi2;
            theta = residualAbsSum(chiNew,sNew,bL,bU,auxdata);
            rho = 0.1;
            nuNew = predPhiDecrease / ((1-rho) * theta);
            prob.nu = max(1,min(1000,nuNew));
                 
            % update predicted and actual reductions with new nu value
            pred = -alphaPr*g*[dchi;ds] - prob.gamma*alphaPr^2*[dchi;ds]'*H*[dchi;ds] + prob.nu*(norm(c',1)-norm(c'+alphaPr*J*[dchi;ds],1));
            ared = merit(prob,chi,chiL,chiU,s,bL,bU,auxdata) - merit(prob,chiNew,chiL,chiU,sNew,bL,bU,auxdata);
        end        

        % Compute acceptance point
        chi = chi + alphaPr * dchi';
        if (ns>=1)
            s = s + alphaPr * ds';
        else
            s = [];
        end
        lambda = lambda + alphaDu * dlambda';
        % Update zLNew and zUNew
        % update from direct solve approach
        zL = zL + alphaDu * dzL';
        zU = zU + alphaDu * dzU';


        % Check for convergence
        sMax = 100; % > 1
        sD = max(sMax,(sum(abs(lambda))+sum(abs(zL))+sum(abs(zU)))/(nc+2*(nchi+ns)));
        sC = max(sMax,(sum(abs(zU))+sum(abs(zL)))/(2*(nchi+ns)));

        part(1) = max(abs(g' + J'*lambda' - zL' + zU'))/sD;
        part(2) = max(abs(lambda'.*c));
        part(3) = max(abs(diag([chi-chiL s-sL])*diag(zL)*e - prob.mu*e))/sC;
        part(4) = max(abs(diag([chiU-chi sU-s])*diag(zU)*e - prob.mu*e))/sC;
        eMu = max(part);        
        % Check for termination conditions
        if (eMu <= prob.eTol)
            fprintf(1,'\nSuccessful solution\n');
            status = 'success';
            break;
        end

        % Check for new barrier problem
        kMu  = 0.8; % (0,1)
        if (eMu < kMu * prob.mu)
            thMu = 1.5; % (1,2)
            % update mu
            prob.mu = max(prob.eTol/10,min(kMu*prob.mu,prob.mu^thMu));
            % update tau
            tau = min(prob.tauMax,100*prob.mu);
        end


        % print iteration
        iprint(prob,iter,chi,lambda,zL,zU,alphaPr,alphaDu,s,chiL,chiU,bL,bU,auxdata);
        

        % reset alpha
        alphaPr = 0.06;
        alphaDu = 0.06;

        % don't do a line search in this MATLAB version
        % just cycle through on another iteration with a lower alpha if
        % the point was not accepted        

    end

    chi = reshape(chi,nStates,N+1)';
    
end