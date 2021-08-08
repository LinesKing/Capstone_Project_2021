% Large-scale Interior Point solution methods
% probOPT Nonlinear Programming Solver sourced from http://www.APMonitor.com and http://www.APOPT.com.

clc;
clear all;
close all;

%  Example Problem
%     min   x2*(5+x1)
%     s.t.  x1 * x2 >= 5
%           x1^2 + x2^2 <=20

%% Basic settings for problem
prob = struct;
prob.nu = 10;  % initial merit function weighting of constraint residuals
prob.mu = 10;  % the initial barrier term

% update mu - used for testing
prob.muUpdate = true;

% Line search criteria
% 	1 = reduction in merit function
% 	2 = simple clipping
%   3 = filter method
prob.lineSearch = 2;

% Solving for z
%   1 = update from direct solve approach
%   2 = update explicitly from z(i) = mu / x(i)
prob.zUpdate = 1;

% reduced or full matrix inversion
% 	1 = condensed, symmetric matrix
% 	2 = full, unsymmetric matrix
prob.matrix = 1;

% show contour plots
prob.contour = true;

% maximum iterations
prob.maxIter = 100;

% show contour plots
prob.showContour = true;

% Initialize slack variables 
% 	true = equation residuals
% 	false = near zero 0.01
prob.slackInit = true;

prob.tauMax = 0.01; % update tau
prob.eTol = 1e-6;  % solution error tolerance
prob.k_mu  = 0.2; % check for new barrier problem (0,1)

%% Primal dual interior method
% Initialize primal variables
x = [2, 3];  % initial value for x
xL = [1, 1];  % lower bound for x
xU = [5, 5];  % upper bound for x
bL = [5, -1e20];  % lower bound for b
bU = [1e20, 20];  % upper bound for b

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
r = res(x,s,bL,bU);

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
% check for consistent bounds
if (min(xU-xL)<0)
    disp('bounds error (xU < xL)')
    sol = [];
    return
end
if (min(bU-bL)<0)
    disp('bounds error (bU < bL)')
    sol = [];
    return
end
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
g = objGrad(x,s);
J = jac(x,bL,bU);

lam = pinv(full(J*J'))*J*(zL'-zU'-g');  % ??
lam = lam';

% Initial parameters
alphaPr = 1.0;  % primal alpha
alphaDu = 1.0;  % dual alpha

% 1-norm of infeasibilities
th = theta(x,s,bL,bU);
thMax = 10^4  * max(1,th);
thMin = 10^-4 * max(1,th);

% Initialize iteration count
iter = 0;

% initialize filter
filter = [thMax phi(prob,x,xL,xU,s,bL,bU)];

% Print iteration zero
iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU);

store = [iter x alphaPr alphaDu prob.mu];

% Start iterating
for iter = 1:prob.maxIter
    
    % new residuals
    r = res(x,s,bL,bU);
    
    % 1-norm of infeasibilities
    th = theta(x,s,bL,bU);
    
    % phi, objective function and barrier terms
    ph = phi(prob,x,xL,xU,s,bL,bU);

    if (prob.zUpdate==2)
        % update explicitly from z = mu / x
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
    end
    
    % make diagonal matrices
    ZL = diag(zL);
    ZU = diag(zU);
    
    % sigmas
    dL = [x-xL s-sL];
    dU = [xU-x sU-s];
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
    J = jac(x,bL,bU);
    W = hes(x,s,lam);
    g = objGrad(x,s);
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
    
    % Compute acceptance point
    xa = x + alphaPr * dx';
    if (ns>=1)
        sa = s + alphaPr * ds';
    else
        sa = [];
    end
    lama = lam + alphaDu * dlam';
    
    % Update zLa and zUa
    switch(prob.zUpdate)
        case(1)
            % update from direct solve approach
            zLa = zL + alphaDu * dzL';
            zUa = zU + alphaDu * dzU';
        case(2)
            % update explicitly from z = mu / x
            for i = 1:n
                zLa(i) = prob.mu / (xa(i)-xL(i));
                dzL(i,1) = zLa(i) - zL(i);
            end
            for i = 1:ns
                zLa(n+i) = prob.mu / (sa(i)-sL(i));
                dzL(n+i,1) = zLa(n+i) - zL(n+i);
            end
            % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
            for i = 1:n
                zUa(i) = prob.mu / (xU(i)-xa(i));
                dzU(i,1) = zUa(i) - zU(i);
            end
            for i = 1:ns
                zUa(n+i) = prob.mu / (sU(i)-sa(i));
                dzU(n+i,1) = zUa(n+i) - zU(n+i);
            end
    end
    
    % max alpha is that which brings the search point to within "tau" of constraint
    % tau is 0 to 0.01 (tau = mu when mu<0.01, otherwise tau=0.01)
    alphaPrMax = 1.0;
    alphaDuMax = 1.0;
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
        if (prob.zUpdate==1)         
            if(zLa(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zL(i)-zL(i))/dzL(i,1));
                %zLa(i) = 0;
            end
            if(zUa(i)<0)
                alphaDuMax = min(alphaDuMax,(tau*zU(i)-zU(i))/dzU(i,1));
                %zUa(i) = 0;
            end
        end
    end
    
    % Line search
    % 	1 = reduction in merit function
    % 	2 = simple clipping
    %   3 = filter method
    switch prob.lineSearch
        case(1)
            % alphaDu set as approach to constraint
            alphaDu = alphaDuMax;
            % compute new acceptance point
            alphaPr = min(alphaPr,alphaPrMax);
        case(2)
            % test case for alpha_pr, alpha_du = 1 with clipping only
            alphaPr = 1;
            alphaDu = 1;
        case(3)
            % alphaDu set as approach to constraint
            alphaDu = alphaDuMax;
            % compute new acceptance point
            alphaPr = min(alphaPr,alphaPrMax);
    end
    
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
    
    % Predicted reduction in the merit function
    pred = -alphaPr*g*[dx;ds] - 0.5*alphaPr^2*[dx;ds]'*H*[dx;ds] + prob.nu*(norm(r',1)-norm(r'+alphaPr*J*[dx;ds],1));
    ared = merit(prob,x,xL,xU,s,bL,bU) - merit(prob,xa,xL,xU,sa,bL,bU);
    
    % Line search criteria
    % 	1 = reduction in merit function
    % 	2 = simple clipping
    % 	3 = filter method
    switch prob.lineSearch
        case(1)
            % merit function
            predPhi1 = g*[dx;ds];
            % set 2nd derivative contribution to zero if < 0
            predPhi2 = max(0,0.5*[dx;ds]'*H*[dx;ds]);
            predPhiDecrease = predPhi1 + predPhi2;
            tha = theta(prob,xa,sa,bL,bU);
            rho = 0.1;
            newNu = predPhiDecrease / ((1-rho) * tha);
            % if (new_nu>prob.nu),
            % 	prob.nu = min(1000,newNu + 1);
            % end
            prob.nu = max(1,min(1000,newNu));
            % update predicted and actual reductions with new nu value
            pred = -alphaPr*g*[dx;ds] - 0.5*alphaPr^2*[dx;ds]'*H*[dx;ds] + prob.nu*(norm(r',1)-norm(r'+alphaPr*J*[dx;ds],1));
            ared = merit(prob,x,xL,xU,s,bL,bU) - merit(prob,xa,xL,xU,sa,bL,bU);
            eta = 0.2;
            % compare actual reduction to predicted reduction
            % as long as the actual reduction is a fraction of the
            % predicted reduction then accept the trial point
            if (ared>=eta*pred)
                ac = true;
            else
                ac = false;
            end
        case(2)
            ac = true;
            % test case for alpha_pr, alpha_du = 1 with clipping only
            alphaPr = 1;
            alphaDu = 1;
        case(3)
            % check if acceptable point with filter method
            [ac] = accept(x,xa,xL,xU);
    end
    
    % optionally show contour plot
    if (prob.contour)
       drawContour(prob,x,alphaPr * dx);
    end
    
    if (ac)
        % accept point
        x = xa;
        s = sa;
        lam = lama;
        zL = zLa;
        zU = zUa;
        
        % update filter
        filter = [filter; th ph];
        %else
        %if (bp_theta(bp,xa)>bp_theta(bp,x)),
        % apply second order correction
        
        %end
    end
    
    % check for convergence
    sMax = 100; % > 1
    sD = max(sMax,(sum(abs(lam))+sum(abs(zL))+sum(abs(zU)))/(m+2*(n+ns)));
    sC = max(sMax,(sum(abs(zU))+sum(abs(zL)))/(2*(n+ns)));
    
    part(1) = max(abs(objGrad(x,s)' + jac(x,bL,bU)'*lam' - zL' + zU'))/sD;
    part(2) = max(abs(res(x,s,bL,bU)));
    % use mu = 0 to test for convergence
    %part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e - bp.mu*e))/sC;
    %part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e - bp.mu*e))/sC;
    part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e))/sC;
    part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e))/sC;
    eMu = max(part);

    % check for termination conditions
    if (eMu <= prob.eTol),
        fprintf(1,'\nSuccessful solution\n');
        status = 'success';
        break;
    end
    
    % Check for new barrier problem
    kMu  = 0.2; % (0,1)

    % Only update mu if requested (not debugging)
    if (prob.muUpdate)
        if (eMu < kMu * prob.mu)
            thMu = 1.5; % (1,2)
            % update mu
            prob.mu = max(prob.eTol/10,min(kMu*prob.mu,prob.mu^thMu));
            % update tau
            tau = min(prob.tauMax,100*prob.mu);
            % re-initialize filter
            filter = [thMax phi(prob,x,xL,xU,s,bL,bU)];
        end
    end

    % print iteration
    iprint(prob,iter,x,lam,zL,zU,alphaPr,alphaDu,s,xL,xU,bL,bU);
    
    store = [store; iter x alphaPr alphaDu prob.mu];

    % reached the end of the iteration loop without convergence
    if (iter==prob.maxIter),
        status = 'failed: max iterations';
        break
    end

    % don't do a line search in this MATLAB version
    % just cycle through on another iteration with a lower alpha if
    % the point was not accepted
    if (ac)
        % reset alpha
        alphaPr = 1.0;
        alphaDu = 1.0;
    else
        % reject point and move alpha_x
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

%% Display figures on iteration progress
figure(2)
hold off;
i = 1;
plot(store(:,1),store(:,i+1),'k-')
hold on;
i = 2;
plot(store(:,1),store(:,i+1),'b-')
legend('x_1','x_2')

figure(3)
hold off;
plot(store(:,1),store(:,n+2),'r-')
hold on;
plot(store(:,1),store(:,n+3),'b-')
plot(store(:,1),log10(store(:,n+4)),'g-')
legend('alpha_{pr}','alpha_{du}','log_{10}(mu)');

% record solution for function return
sol = struct;
sol.status = status;
sol.x = x;
sol.s = s;
sol.lam = lam;
sol.zL = zL;
sol.zU = zU;