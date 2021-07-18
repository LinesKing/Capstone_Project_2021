function [dxidt] = unicycleODE(t, xi, u)
%%%% This function is for unicycle model, modified so angular velocity and
%%%% acceleration are the inputs
%%%% Input: t, xi, u
%%%% t: times (required for ode solver)
%%%% xi: states of system [x, y, v, theta]
%%%% u: inputs in system [omega, a]
%%%% Output: delta_u, delta_xi
%%%%    dxidt: time derivative of states

    dxidt = zeros(4,1);
    dxidt(1) = xi(4)*cos(xi(3));
    dxidt(2) = xi(4)*sin(xi(3));
    dxidt(3) = u(1);
    dxidt(4) = u(2);
    
end

