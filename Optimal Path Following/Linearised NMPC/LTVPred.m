function xiNext = LTVPred(xi, u, T)
%%%% This function is to predict next state by present state using a LTV
%%%% model.
%%%% Input: xi, u ,T
%%%%    xi: Current states.
%%%%    u: Current inputs.
%%%%    T: Sampled time.
%%%% Output: xiNext
%%%%    xiNext: Next states.

    dx = xi(4)*cos(xi(3));
    ddx = -xi(4)*sin(xi(3))*u(1) + cos(xi(3))*u(2);
    dy = xi(4)*sin(xi(3));
    ddy = xi(4)*cos(xi(3))*u(1) + sin(xi(3))*u(2);
    dtheta = u(1);
    dv = u(2);
    
    xiNext = zeros(4,1);
    xiNext(1) = xi(1) + T*dx + T^2/2*ddx;
    xiNext(2) = xi(2) + T*dy + T^2/2*ddy;
    xiNext(3) = xi(3) + T*dtheta;
    xiNext(4) = xi(4) + T*dv;
end

