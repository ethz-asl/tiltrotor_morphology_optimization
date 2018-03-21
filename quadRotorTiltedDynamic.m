function [pdotdot, wbdot] = quadRotorTiltedDynamic(kf, km, Ib, wRb, alpha,n, L, g, m)
%function [pdotdot, wbdot] = quadRotorTiltedDynamic(kf, km, Ib, wRb, alpha,n)
%QUADROTORTILTEDDYNAMIC returns the dynamic of a quadcopter with tilted
%propellers. Returns the linear and angular acceleration of the drone.
f = [0 0 -g].';
F = [0 kf*sin(alpha(2)) 0 -kf*sin(alpha(4)); -kf*sin(alpha(1)) 0 kf*sin(alpha(3)) 0;
    kf*cos(alpha(1)) kf*cos(alpha(2)) kf*cos(alpha(3)) kf*cos(alpha(4))];

Tau = [0 L*kf*cos(alpha(2))-km*sin(alpha(2)) 0 -L*kf*cos(alpha(4))+km*sin(alpha(4));
    -L*kf*cos(alpha(1))-km*sin(alpha(1)) 0 L*kf*cos(alpha(3))+km*sin(alpha(3)) 0;
    -L*kf*sin(alpha(1))+km*cos(alpha(1)) -L*kf*sin(alpha(2))-km*cos(alpha(2)) -L*kf*sin(alpha(3))+km*cos(alpha(3)) -L*kf*sin(alpha(4))-km*cos(alpha(4))];


pdotdot = f + F*n.^2/m;

wbdot = Ib\(Tau*n.^2);

Ndecimals = 8;
k = 10.^Ndecimals;
pdotdot = round(k*pdotdot)/k;
wbdot = round(k*wbdot)/k;

end

