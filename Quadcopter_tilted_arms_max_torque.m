function [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, nhover, alpha0, n0, d, beta, theta)
% [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, nmax, nhover, alpha0, n0, d, beta, theta)
%QUADCOPTER_TILTED_ARMS_MAX_TORQUE Summary of this function goes here
%   Detailed explanation goes here

%% Optimization of alpha and n
% maximize norm of the Torque M in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% Condition  Ax <= b        
A = []; 
b = [];                 

% condition: Aeq.x = beq     
Aeq = [];
beq = [];

% condition: lb <= x <= ub 
lb = [-pi -pi -pi -pi 0 0 0 0].'; % lower bound
ub = [pi pi pi pi nmax nmax nmax nmax].'; % upper bound

% initial guess:
x0 = [alpha0, n0.'];

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
x = fmincon(@(x) TorqueNorm(kf, km, beta, theta, L, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconM(x, kf, km, beta, theta, L, d),  options);

% Substitution:
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
k = 10.^Ndecimals;
nstar = round(k*nstar)/k;
alphastar = round(k*alphastar)/k;


end
function [squareNormM] = TorqueNorm(kf, km, beta, theta, L, x)
M = Torque(kf, km, theta, beta, L, x);
% Objecive function maximize the squared norm of the torque M:
squareNormM = -(M(1)^2+ M(2)^2+ M(3)^2);
end
function [M] = Torque(kf, km, theta, beta, L, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
M1 = km*x(5)^2*(sin(x(1))*sin(theta(1)) + cos(x(1))*sin(beta(1))*cos(theta(1))) - km*x(8)^2*(sin(x(4))*sin(theta(4) + (3*pi)/2) + cos(x(4))*sin(beta(4))*cos(theta(4) + (3*pi)/2)) - km*x(6)^2*(sin(x(2))*sin(theta(2) + pi/2) + cos(x(2))*sin(beta(2))*cos(theta(2) + pi/2)) - km*x(7)^2*(sin(x(3))*sin(theta(3)) + cos(x(3))*sin(beta(3))*cos(theta(3))) - L*kf*x(5)^2*sin(beta(1))*sin(x(1))*cos(theta(1)) + L*kf*x(7)^2*sin(beta(3))*sin(x(3))*cos(theta(3)) - L*kf*x(6)^2*sin(beta(2))*sin(x(2))*cos(theta(2) + pi/2) - L*kf*x(8)^2*sin(beta(4))*sin(x(4))*cos(theta(4) + (3*pi)/2);
M2 = km*x(6)^2*(sin(x(2))*cos(theta(2) + pi/2) - cos(x(2))*sin(beta(2))*sin(theta(2) + pi/2)) + km*x(8)^2*(sin(x(4))*cos(theta(4) + (3*pi)/2) - cos(x(4))*sin(beta(4))*sin(theta(4) + (3*pi)/2)) - km*x(5)^2*(sin(x(1))*cos(theta(1)) - cos(x(1))*sin(beta(1))*sin(theta(1))) + km*x(7)^2*(sin(x(3))*cos(theta(3)) - cos(x(3))*sin(beta(3))*sin(theta(3))) - L*kf*x(5)^2*sin(beta(1))*sin(x(1))*sin(theta(1))+ L*kf*x(7)^2*sin(beta(3))*sin(x(3))*sin(theta(3)) - L*kf*x(6)^2*sin(beta(2))*sin(x(2))*sin(theta(2) + pi/2) - L*kf*x(8)^2*sin(beta(4))*sin(x(4))*sin(theta(4) + (3*pi)/2);
M3 = km*x(5)^2*cos(x(1))*cos(beta(1)) - km*x(6)^2*cos(x(2))*cos(beta(2)) + km*x(7)^2*cos(x(3))*cos(beta(3)) - km*x(8)^2*cos(x(4))*cos(beta(4)) - L*kf*x(5)^2*cos(beta(1))*sin(x(1)) -L*kf*x(6)^2*cos(beta(2))*sin(x(2)) - L*kf*x(7)^2*cos(beta(3))*sin(x(3)) - L*kf*x(8)^2*cos(beta(4))*sin(x(4));
M = [M1 M2 M3];
end
function [c,ceq] = nonlinconM(x, kf, km, beta, theta,  L, d)
% function [c,ceq] = nonlincon(x, kf, nmax,d, L)

c = [];

%Torque
M = Torque(kf, km, theta, beta, L, x);

% Torque vector parallel to d condition:
ceq(1) = M(2)*d(3) - M(3)*d(2);
ceq(2) = M(3)*d(1) - M(1)*d(3);
ceq(3) = M(1)*d(2) - M(2)*d(1);
% % Force applied on MAV equal to zero condition:
% ceq(4) = T(1);
% ceq(5) = T(2);
% ceq(6) = T(3);
end

