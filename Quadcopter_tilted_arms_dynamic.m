function [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, teta,n, L, g, Mb, Mp, R)
%function [pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, teta,n, L, g, Mb, Mp, R)
%QUADROTORTILTEDDYNAMIC returns the dynamic of a quadcopter with tilting
%propeller and tilted arms
%propellers. Returns the linear and angular acceleration of the drone, its inertia tensor and mass.
bRp1 = rotz(rad2deg(teta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
bRp2 = rotz(rad2deg(pi/2+teta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
bRp3 = rotz(rad2deg(pi+teta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
bRp4 = rotz(rad2deg(3*pi/2+teta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
Op1 = rotz(rad2deg(teta(1)))*roty(rad2deg(beta(1)))*[L 0 0].'
Op2 = rotz(rad2deg(pi/2+teta(2)))*roty(rad2deg(beta(2)))*[L 0 0].'
Op3 = rotz(rad2deg(pi+teta(3)))*roty(rad2deg(beta(3)))*[L 0 0].'
Op4 = rotz(rad2deg(3*pi/2+teta(4)))*roty(rad2deg(beta(4)))*[L 0 0].'

m = Mb + 4*Mp; % Mass total of the drone
Icom = 2/5*Mb*R*R*eye(3); % Inertia tensor of a sphere
% Inertia tensor of a sphere with rotors represented as point masses
Ip = Mp*(norm(Op1)^2*eye(3) - Op1*Op1.' + norm(Op2)^2*eye(3) - Op2*Op2.' + ...
     norm(Op3)^2*eye(3) - Op3*Op3.' + norm(Op4)^2*eye(3) - Op4*Op4.');% Inertia tensor of rotors (point masses)
Ib = Icom+ Ip; % Inertia tensor of the drone (sphere with 4 point masses)
 
f = [0 0 -g].'; % gravity

Tp1 = [0 0 kf*n(1)^2].'; % Thrust vector propeller 1
Tp2 = [0 0 kf*n(2)^2].'; % Thrust vector propeller 2
Tp3 = [0 0 kf*n(3)^2].'; % Thrust vector propeller 3
Tp4 = [0 0 kf*n(4)^2].'; % Thrust vector propeller 4
Tauext1 = [0 0 -km*n(1)^2].'; % Thrust vector propeller 1
Tauext2 = [0 0 km*n(2)^2].'; % Thrust vector propeller 2
Tauext3 = [0 0 -km*n(3)^2].'; % Thrust vector propeller 3
Tauext4 = [0 0 km*n(4)^2].'; % Thrust vector propeller 4
Taub = cross(Op1,bRp1*Tp1) + cross(Op2,bRp2*Tp2) + cross(Op3,bRp3*Tp3) + cross(Op4,bRp4*Tp4);
pdotdot = f + (1/m)*wRb*(bRp1*Tp1 + bRp2*Tp2 + bRp3*Tp3 + bRp4*Tp4);

wbdot = Ib\(bRp1*Tauext1 + bRp2*Tauext2 + bRp3*Tauext3 + bRp4*Tauext4 + Taub);

Ndecimals = 8;
k = 10.^Ndecimals;
pdotdot = round(k*pdotdot)/k;
wbdot = round(k*wbdot)/k;

end