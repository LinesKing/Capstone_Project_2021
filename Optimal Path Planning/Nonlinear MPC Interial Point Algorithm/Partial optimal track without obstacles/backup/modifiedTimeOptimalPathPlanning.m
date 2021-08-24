function [x] = modifiedTimeOptimalPathPlanning(xInner, yInner, xOuter, yOuter, warm, initPoint)
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
    prob.mu = 1;  % the initial barrier term

    % update mu 
    prob.muUpdate = true;

    % maximum iterations
    prob.maxIter = 1000;

    % Initialize slack variables 
    % 	true = equation residuals
    % 	false = near zero 0.01
    prob.slackInit = true;

    prob.tauMax = 0.01; % update tau
    prob.eTol = 1e-1;  % solution error tolerance
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
    weightDist = 0.1;
    weightVel = 0.5; % 0.5
    weightTime = 1; % 1

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
        x = [x, warm.xWarm(i), warm.yWarm(i), warm.vxWarm(i), warm.vyWarm(i), warm.axWarm(i), warm.ayWarm(i), warm.tWarm];

        % Upper and lower bounds of decision variables (box constraints)
        % [x, y, vx, vy, ax, ay]
        if (i == 1) % satisfy initial conditions contraint
            xL = [xL, [initPoint(1) initPoint(2) vxMin vyMin axMin ayMin]];
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
        xL = [xL, tMin];
        xU = [xU, tMax];
    end
    xL = xL - 1e-5;
    xU = xU + 1e-5;
    
    % Initialise upper and lower bounds for constraints
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

    tau = min(prob.tauMax,100*prob.mu);  % ??

    % Variable constraint multipliers
    % zL*(x-xL) = mu  =>  zL = mu / (x-xL)
    for i = 1:n
        zL(i) = prob.mu / (x(i)-xL(i));
    end
    for i = 1:ns
        zL(n+i) = prob.mu / (s(i)-sL(i));
    end
    % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
    for i = 1:n
        zU(i) = prob.mu / (xU(i)-x(i));
    end
    for i = 1:ns
        zU(n+i) = prob.mu / (sU(i)-s(i));
    end

    % Initialize equation constraint multipliers 
    g = objGrad(x,s,auxdata);
    J = jac(x,bL,bU,auxdata);

%     lam = pinv(full(J*J'))*J*(zL'-zU'-g'); 
    lam = lsqminnorm(full(J*J'),J)*(zL'-zU'-g');
    lam = lam';

    % Initial parameters
    alphaPr = 0.02;  % primal alpha
    alphaDu = 0.02;  % dual alpha

    % Initialize iteration count
    iter = 0;

    % Print iteration zero
    iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU,auxdata);

    % Start iterating
    for iter = 1:prob.maxIter

        % new residuals
        r = res(x,s,bL,bU,auxdata);

        % phi, objective function and barrier terms
        ph = phi(prob,x,xL,xU,s,bL,bU,auxdata);

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
        % d = -pinv(full(A))*b;
        d = lsqminnorm(A,-b);
        dx = d(1:n);
        ds = d(n+1:n+ns,1);
        dlam = d(n+ns+1:n+ns+m,1);
        % compute search direction for z (explicit solution)
        dzL = prob.mu*invDL*e - zL' - SigL*[dx; ds];
        dzU = prob.mu*invDU*e - zU' + SigU*[dx; ds];
        search = [d;dzL;dzU];
    
        % Compute acceptance point
        xa = x + alphaPr * dx';
        if (ns>=1)
            sa = s + alphaPr * ds';
        else
            sa = [];
        end
        lama = lam + alphaDu * dlam';

        % Update zLa and zUa
        % update from direct solve approach
        zLa = zL + alphaDu * dzL';
        zUa = zU + alphaDu * dzU';

        % max alpha is that which brings the search point to within "tau" of constraint
        % tau is 0 to 0.01 (tau = mu when mu<0.01, otherwise tau=0.01)
        alphaPrMax = 0.02;
        alphaDuMax = 0.02;
        % Check for constraint violations
        for i = 1:n
            if(xa(i)<xL(i))
                alphaPrMax = min(alphaPrMax,(xL(i)+tau*(x(i)-xL(i))-x(i))/dx(i,1));
            end
            if(xa(i)>xU(i))
                alphaPrMax = min(alphaPrMax,(xU(i)-tau*(xU(i)-x(i))-x(i))/dx(i,1));
            end
        end
        for i = 1:ns
            if(sa(i)<sL(i))
                alphaPrMax = min(alphaPrMax,(sL(i)+tau*(s(i)-sL(i))-s(i))/ds(i,1));
            end
            if(sa(i)>sU(i))
                alphaPrMax = min(alphaPrMax,(sU(i)+tau*(sU(i)-s(i))-s(i))/ds(i,1));
            end
        end
        for i = 1:n       
            if(zLa(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zL(i)-zL(i))/dzL(i,1));
                zLa(i) = 0;
            end
            if(zUa(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zU(i)-zU(i))/dzU(i,1));
                zUa(i) = 0;
            end
        end

        % Line search
%         % test case for alpha_pr, alpha_du = 1 with clipping only
%         % alphaDu set as approach to constraint
%         alphaDu = alphaDuMax;
%         % compute new acceptance point
%         alphaPr = min(alphaPr,alphaPrMax);
%         
%         % Predicted reduction in the merit function
%         pred = -alphaPr*g*[dx;ds] - 0.5*alphaPr^2*[dx;ds]'*H*[dx;ds] + prob.nu*(norm(r',1)-norm(r'+alphaPr*J*[dx;ds],1));
%         ared = merit(prob,x,xL,xU,s,bL,bU,auxdata) - merit(prob,xa,xL,xU,sa,bL,bU,auxdata);  
%         
%         % merit function
%         predPhi1 = g*[dx;ds];
%         % set 2nd derivative contribution to zero if < 0
%         predPhi2 = max(0,0.5*[dx;ds]'*H*[dx;ds]);
%         predPhiDecrease = predPhi1 + predPhi2;
%         tha = theta(xa,sa,bL,bU,auxdata);
%         rho = 0.1;
%         newNu = predPhiDecrease / ((1-rho) * tha);
%         % if (new_nu>prob.nu),
%         % 	prob.nu = min(1000,newNu + 1);
%         % end
%         prob.nu = max(1,min(1000,newNu));
%         % update predicted and actual reductions with new nu value
%         pred = -alphaPr*g*[dx;ds] - 0.5*alphaPr^2*[dx;ds]'*H*[dx;ds] + prob.nu*(norm(r',1)-norm(r'+alphaPr*J*[dx;ds],1));
%         ared = merit(prob,x,xL,xU,s,bL,bU,auxdata) - merit(prob,xa,xL,xU,sa,bL,bU,auxdata);
%         eta = 0.2;
%         % compare actual reduction to predicted reduction
%         % as long as the actual reduction is a fraction of the
%         % predicted reduction then accept the trial point
%         if (ared>=eta*pred)
%             ac = true;
%         else
%             ac = false;
%         end
        
        
        c1 = 1e-3; c2 = 0.8;
        while 1
            Armijo = (obj(x+alphaPr*[dx]',auxdata) <= obj(x,auxdata) + c1*alphaPr*g*[dx;ds]);
            % curvature = (objGrad(x+alphaPr*dx,s+alphaPr*ds,auxdata)*[dx;ds] >= c2*g*[dx;ds]);
            if (Armijo)
                ac = true;
            else
                ac = false;
            end
            
            if (ac)
                break
            else
                % reject point and move alpha_x
                alphaPr = alphaPr * 0.5;
                alphaDu = alphaDu * 0.5;
    %             if (alphaPr < 1e-6)
    %                 alphaPr = 0.02;
    %                 alphaDu = 0.02;
    %             end
            end
        end
                
        

        
%         ac = true;

        % Clipping (this should already be arranged by alpha_max)
        % push away from the boundary with tau
        for i = 1:n
            if(xa(i)<xL(i))
                xa(i) = xL(i)+tau*(x(i)-xL(i));
            end
            if(xa(i)>xU(i))
                xa(i) = xU(i)-tau*(xU(i)-x(i));
            end
            if(zLa(i)<0)
                zLa(i) = tau*(zL(i));
            end
            if(zUa(i)<0)
                zUa(i) = tau*(zU(i));
            end
        end
        for i = 1:ns
            if(sa(i)<sL(i))
                sa(i) = sL(i)+tau*(s(i)-sL(i));
            end
            if(sa(i)>sU(i))
                sa(i) = sU(i)-tau*(sU(i)-s(i));
            end
        end


        if (ac)
            % accept point
            x = xa;
            s = sa;
            lam = lama;
            zL = zLa;
            zU = zUa;
        end

        % check for convergence
%         sMax = 10000; % > 1
%         sD = max(sMax,(sum(abs(lam))+sum(abs(zL))+sum(abs(zU)))/(m+2*(n+ns)));
%         sC = max(sMax,(sum(abs(zU))+sum(abs(zL)))/(2*(n+ns)));
% 
%         part(1) = max(abs(objGrad(x,s,auxdata)' + jac(x,bL,bU,auxdata)'*lam' - zL' + zU'))/sD;
%         part(2) = max(abs(res(x,s,bL,bU,auxdata)))/sD*100;
%         % use mu = 0 to test for convergence
%         part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e - prob.mu*e))/sC;
%         part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e - prob.mu*e))/sC;
% %         part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e))/sC;
% %         part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e))/sC;
%         eMu = max(part);
        
        eMu = max(abs(objGrad(x,s,auxdata)' + jac(x,bL,bU,auxdata)'*lam' - zL' + zU'));

%         fprintf(1,'  %11.4e %11.4e %11.4e %11.4e \n',part(1),part(2),part(3),part(4));
        
        % check for termination conditions
        if (eMu <= prob.eTol)
            fprintf(1,'\nSuccessful solution\n');
            status = 'success';
            break;
        end

        % Check for new barrier problem
        kMu  = 0.2; % (0,1)

        % Only update mu if requested
%         if (prob.muUpdate)
%             if (eMu < kMu * prob.mu)
%                 thMu = 1.5; % (1,2)
%                 % update mu
%                 prob.mu = max(prob.eTol/10,min(kMu*prob.mu,prob.mu^thMu));
%                 % update tau
%                 tau = min(prob.tauMax,100*prob.mu);
%             end
%         end
        prob.mu = 0.95*prob.mu;

        % print iteration
        iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU,auxdata);
        
        if (ac)
            % reset alpha
            alphaPr = 0.02;
            alphaDu = 0.02;
        end
        % don't do a line search in this MATLAB version
        % just cycle through on another iteration with a lower alpha if
        % the point was not accepted
        

    end

    x = reshape(x,nStates,N+1)';
    
end