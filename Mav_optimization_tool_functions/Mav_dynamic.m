function [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(n, kf, km, wRb, alpha, beta, theta,w, L, g, dec, gravitiy)
%MAV_DYNAMIC returns the dynamic of a mav with tilting
%propeller and tilted arms
%   Returns the linear and angular acceleration of the drone, its inertia tensor and mass.

interval = 2*pi/n; % interval between arms in normal n-copter configuration
% pre-allocation:
bRp = zeros(3,3,n);
Op = zeros(3,n);
Tp = zeros(3,n);
Tauext = zeros(3,n);
Taub = [0; 0; 0];
M = [0; 0; 0];
T = [0; 0; 0];
for i =1:n
    %% Find the propellers rotation matrix
    bRp(:,:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i))*Rotx(alpha(i));
    
    %% Find the propellers positions in the body frame
    Op(:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i))*[L 0 0].';
    
    %% Find forces applied by all propellers thrusts on the body
    Tp(:,i) = [0 0 kf*w(i)^2].'; % force applied by every propeller in propeller frame
    T  = T + bRp(:,:,i)*Tp(:,i); % force applied by all the propellers in body frame
    
    %% Find torques applied by all propellers on the body (in propeller frame)
    if mod(i,2) == 0
        c = 1; % counter torque defined negative for clock wise rotating propellers (all odd propellers)
    else
        c = -1;
    end
    Tauext(:,i) = [0 0 c*km*w(i)^2].'; % counter torque produced by every propeller in propeller frame
    M = M + bRp(:,:,i)*Tauext(:,i); % torque applied by every propeller on the body
    Taub = Taub + cross(Op(:,i),bRp(:,:,i)*Tp(:,i)); % total torque applied on the body (in body frame)
end
%% Drone inertia
[m, Ib] = Mav_inertias(n, L, theta, beta);

%% take gravity into account?
if gravitiy
    f = [0 0 -g].'; % gravity
else
    f = [0 0 0].';
end
%% angular acceleration in the body frame
M = M + Taub;
wbdot = Ib\M;
wbdot = round(dec*wbdot)/dec; % rounded

%% linear acceleration in the body frame
pdotdot = f + (1/m)*wRb*T;
pdotdot = round(dec*pdotdot)/dec; %rounded
end