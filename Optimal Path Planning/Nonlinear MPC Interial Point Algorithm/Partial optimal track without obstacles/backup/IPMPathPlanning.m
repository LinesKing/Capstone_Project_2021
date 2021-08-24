function [x] = IPMPathPlanning(xInner, yInner, xOuter, yOuter, warm, initPoint)
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
    prob.nu = 10;  % initial merit function weighting of constraint residuals
    prob.mu = 10;  % the initial barrier term

    % maximum iterations
    prob.maxIter = 500;

    % Initialize slack variables 
    % 	true = equation residuals
    % 	false = near zero 0.01
    prob.slackInit = true;

    prob.tauMax = 0.01; % update tau
    prob.eTol = 1e-7;  % solution error tolerance

    
    % Receeding horizon
    N = length(xOuter) - 1;

    % Max time between points (Must be chosen to something physcally poissble)
    tMax = 0.05; 

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
    weightDist = 0.1;
    weightVel = 0.5;
    weightTime = 1;

    % Define constants (number of variables)
    nSystemStates = 4;
    nControlInputs = 2;
    nTimeVariable = 1;
    nStates = nSystemStates + nControlInputs + nTimeVariable;

    % Set up the auxiliary data
    auxdata = {N nStates weightDist weightVel weightTime ...
                   xOuter yOuter xInner yInner};
    %% Primal dual interior method
    % Initialise decision variables, lower bounds and upper bounds
    x = [];
    xL = [];
    xU = [];
    for i = 1:N+1
        % Set initial start point to center of track (warm start)
        % [x, y, vx, vy, ax, ay, t]
        x = [x, warm.xWarm(i), warm.yWarm(i), zeros(1, nStates-2)];

        % Upper and lower bounds of decision variables (not constraint equations)
        % [x, y, vx, vy, ax, ay]
        if (i == 1) % satisfy initial conditions contraint
            xL = [xL, [initPoint(1) initPoint(2) vxMin vyMin axMin ayMin]]; %2D array with each row filled with column of vector current states
            xU = [xU, [initPoint(1) initPoint(2) vxMax vyMax axMax ayMax]];  
        else
            yLb = min(yOuter(i),yInner(i));
            yUb = max(yOuter(i),yInner(i));
            xLb = min(xOuter(i),xInner(i));
            xUb = max(xOuter(i),xInner(i));
            xL = [xL, [xLb yLb vxMin vyMin axMin ayMin]];
            xU = [xU, [xUb yUb vxMax vyMax axMax ayMax]]; 
        end

        % Lower and Upper Bounds for sampling time
        xL = [xL, 0];
        xU = [xU, tMax];
    end

    % Initialise upper and lower bounds for constraints
    bL = [];
    bU = [];
    for i = 1:N
        % System Constraints
        bL = [bL, [0 0 0 0]];
        bU = [bU, [0 0 0 0]];

        % Triplet contraints
        bL = [bL, 0];
        bU = [bU, 0];
    end    

    % Add slack for inequality constraints only
    si = [];  % iter for inequality constraints
    s = [];  % slack variables
    sL = [];  % lower bound for s
    sU = [];  % upper bound for s
    k = 0;
    for i = 1:max(size(bL))
        if(bU(i)>bL(i))
            % inequality constraints
            k = k + 1;
            si(k) = i;
            s(k) = 0;
            sL(k) = bL(i);
            sU(k) = bU(i);
        end
    end

    % Problem size
    n = max(size(x));  % number of primal variables x
    ns = max(size(s));  % number of slack variables s

    % Initialize residuals
    r = res(x,s,bL,bU,auxdata);

    % number of constraints m (inequality or equality)
    m = max(size(r));
    % ones vector
    e = ones(n+ns,1);

    % Initial equation slack variables
    k = 0;
    for i = 1:m
        % no slack variable when bU(i) == bL(i)
        if(bU(i)>bL(i))
            k = k + 1;
            if (prob.slackInit)
               s(k) = r(i);
            else
               s(k) = 0.01;
            end
        end
    end

    % Check and make the problem feasible
    % move x into feasible region and off boundary initially
    for i = 1:n
        if (x(i)<=xL(i) || x(i)>=xU(i))
            fprintf(1,'Moving x0 to interior region\n')
            break
        end
    end
    % move x variables to be feasible
    for i = 1:n
        if (x(i)<=xL(i))
            x(i) = min(xU(i),xL(i)+1e-2);
        end
        if (x(i)>=xU(i))
            x(i) = max(xL(i),xU(i)-1e-2);
        end
    end
    % move slack variables to be feasible
    for i = 1:ns
        if (s(i)<=sL(i))
            s(i) = min(sU(i),sL(i)+1e-2);
        end
        if (s(i)>=sU(i))
            s(i) = max(sL(i),sU(i)-1e-2);
        end
    end

    % Variable constraint multipliers
    % zL*(x-xL) = mu  =>  zL = mu / (x-xL)
    for i = 1:n
        zL(i) = min(prob.mu / (x(i)-xL(i)),10e8);
    end
    for i = 1:ns
        zL(n+i) = min(prob.mu / (s(i)-sL(i)),10e8);
    end
    % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
    for i = 1:n
        zU(i) = min(prob.mu / (xU(i)-x(i)),10e8);
    end
    for i = 1:ns
        zU(n+i) = min(prob.mu / (sU(i)-s(i)),10e8);
    end

    % Initialize equation constraint multipliers 
    g = objGrad(x,s,auxdata);
    J = jac(x,bL,bU,auxdata);

    lam = pinv(full(J*J'))*J*(zL'-zU'-g');  % ??
    lam = lam';

    % Initial parameters
    alphaPr = 1;  % primal alphaPr
    alphaDu = 1;  % dual alphaPr

    % Initialize iteration count
    iter = 0;

    % Print iteration zero
    iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU,auxdata);

    store = [iter x alphaPr alphaDu prob.mu];

    % Start iterating
    for iter = 1:prob.maxIter

        % new residuals
        r = res(x,s,bL,bU,auxdata);

        % make diagonal matrices
        ZL = diag(zL);
        ZU = diag(zU);

        % sigmas
        dL = [x-xL+10-8 s-sL+10-8];
        dU = [xU-x+10-8 sU-s+10-8];
        DL = diag(dL);
        DU = diag(dU);
        invDL = zeros(n+ns,n+ns);
        invDU = zeros(n+ns,n+ns);
        SigL = zeros(n+ns,n+ns);
        SigU = zeros(n+ns,n+ns);
        for i = 1:n+ns
           invDL(i,i) = 1 / DL(i,i); 
           invDU(i,i) = 1 / DU(i,i); 
           SigL(i,i) = ZL(i,i) / DL(i,i);
           SigU(i,i) = ZU(i,i) / DU(i,i);
        end

        % 1st and 2nd derivatives
        J = jac(x,bL,bU,auxdata);
        W = hes(x,s,lam,auxdata);
        g = objGrad(x,s,auxdata);
        H = W + SigL + SigU;

        % Construct and solve Ax=b
        % Symmetric, zL/zU explicitly solved  ??
        % construct A and b
        A = [H,J';J,zeros(m,m)];
        % b(1:n+ns,1) = g' - zL' + zU' + J'*lam';
        b(1:n+ns,1) = g' + J'*lam' - prob.mu*invDL*e + prob.mu*invDU*e;
        b(n+ns+1:n+ns+m,1) = r(1:m)';
        % compute search direction, solving Ax = b
        d = -pinv(full(A))*b;
        dx = d(1:n,1);
        ds = d(n+1:n+ns,1);
        dlam = d(n+ns+1:n+ns+m,1);
        % compute search direction for z (explicit solution)
        dzL = prob.mu*invDL*e - zL' - SigL*[dx; ds];
        dzU = prob.mu*invDU*e - zU' + SigU*[dx; ds];
        search = [d;dzL;dzU];
        
        % Line search
        gamma = 1e-4;
        delta = 0.5;
        rhok  = 1e-8;
        
        % Primal line search
        while (phi(bp,x+alphaPr.*dx,xL,xU,s+alphaPr.*ds,bL,bU,auxdata)>phi(bp,x,xL,xU,s,bL,bU,auxdata)-gamma*alphaPr^2*(norm([dx ds]))^2) 
            if (alphaPr*norm([dx ds]) < rhok)   
                alphaPr  = 0;              % <-- failure to search for a value of alphaPr nonzero
            else
                alphaPr = alphaPr*delta;     % <-- reduction of the steplength
            end
        end 

        
        % Dual line search
        while (theta(x+alphaDu.*dk,s,bL,bU,auxdata)>theta(x,s,bL,bU,auxdata)-gamma*alphaDu^2*(norm(dk))^2) 
            if (alphaDu*norm(dk) < rhok)   
                alphaDu  = 0;              % <-- failure to search for a value of alpha nonzero
            else
                alphaDu = alphaDu*delta;     % <-- reduction of the steplength
            end
        end 
        alpha1 = alphaDu;
        F1     = theta(x+alpha1.*dk,s,bL,bU,auxdata)-(theta(x,s,bL,bU,auxdata)-gamma*alpha1^2*(norm(dk))^2);

        
        % Compute acceptance point
        x = x + alphaPr * dx';
        if (ns>=1)
            s = s + alphaPr * ds';
        else
            s = [];
        end
        lam = lam + alphaDu * dlam';

        % Update zLa and zUa
        % update from direct solve approach
        zL = zL + alphaDu * dzL';
        zU = zU + alphaDu * dzU';

        % check for convergence
        sMax = 100; % > 1
        sD = max(sMax,(sum(abs(lam))+sum(abs(zL))+sum(abs(zU)))/(m+2*(n+ns)));
        sC = max(sMax,(sum(abs(zU))+sum(abs(zL)))/(2*(n+ns)));

        part(1) = max(abs(objGrad(x,s,auxdata)' + jac(x,bL,bU,auxdata)'*lam' - zL' + zU'))/sD;
        part(2) = max(abs(res(x,s,bL,bU,auxdata)));
        % use mu = 0 to test for convergence
        %part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e - bp.mu*e))/sC;
        %part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e - bp.mu*e))/sC;
        part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e))/sC;
        part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e))/sC;
        eMu = max(part);

        % check for termination conditions
        if (eMu <= prob.eTol)
            fprintf(1,'\nSuccessful solution\n');
            status = 'success';
            break;
        end

        % Check for new barrier problem
        kMu  = 0.2; % (0,1)

        % Only update mu if requested (not debugging)
        if (eMu < kMu * prob.mu)
            thMu = 1.5; % (1,2)
            % update mu
            prob.mu = max(prob.eTol/10,min(kMu*prob.mu,prob.mu^thMu));
        end


        % print iteration
        iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU,auxdata);

        store = [store; iter x alphaPr alphaDu prob.mu];

        % reached the end of the iteration loop without convergence
        if (iter==prob.maxIter)
            status = 'failed: max iterations';
            break
        end

        % don't do a line search in this MATLAB version
        % just cycle through on another iteration with a lower alphaPr if
        % the point was not accepted
        if (ac)
            % reset alphaPr
            alphaPr = 1.0;
            alphaDu = 1.0;
        else
            % reject point and move alphaPr_x
            alphaPr = alphaPr / 2;
            alphaDu = alphaDu / 2;
            if (alphaPr < 1e-4)
                alphaPr = 1.0;
                alphaDu = 1.0;
            end
        end
    end

    %% Display final solution
    disp('Status')
    disp(status)
    disp('Solution: ')
    disp(x)
    disp('Slack Variables: ')
    disp(s)
    disp('Equation Multipliers: ')
    disp(lam)
    disp('Lower Constraint Variable Mult: ')
    disp(zL)
    disp('Upper Constraint Variable Mult: ')
    disp(zU)

    %% Record solution for function return
    sol = struct;
    sol.status = status;
    sol.x = x;
    sol.s = s;
    sol.lam = lam;
    sol.zL = zL;
    sol.zU = zU;

end