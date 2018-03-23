%% Parameters
M = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
m = M + 4*Mp; % total mass [kg]
g = 9.81;
kf = 3.86e-4; % Propeller thrust coefficient
km = 4e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
%% Inertia
Op1 = [L 0 0].';
Op2 = [0 L 0].';
Op3 = [-L 0 0].';
Op4 = [0 -L 0].';
Icom = 2/5*m*R*R*eye(3);
% Inertia tensor of a spher with 
Ib = Mp*(norm(Op1)^2*eye(3) - Op1*Op1.' + norm(Op2)^2*eye(3) - Op2*Op2.' + ...
     norm(Op3)^2*eye(3) - Op3*Op3.' + norm(Op4)^2*eye(3) - Op4*Op4.') +Icom;
%% Init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb = rotz(roll0)*roty(pitch0)*rotz(yaw0);
alpha = [0 0 0 0].';
n = [95 95 96 96].';
wb = [0 0 0].';
pdot = [0 0 0].';
p = [0 0 1].';
h = p;
%% Simulation

% Simulation time [s]
start_time = 0;
end_time = 2;
dt = 0.005;
times = start_time:dt:end_time;

% Step through the simulation, updating the state.
for t = times

    % Compute linear and angular accelerations.
    [a, wbdot] = quadRotorDynamic(kf, km, Ib, wRb, alpha,n, L,g, m);
    wb = wb + dt * wbdot;
    wb = mod(wb,2*pi);
    wRbdot = wRb*skew(wb);
    wRb = wRb + dt*wRbdot;
    pdot = pdot + dt * a;
    p = p + dt * pdot;
    h = [h p];
end
figure (1);
plot3(h(1,:),h(2,:),h(3,:));