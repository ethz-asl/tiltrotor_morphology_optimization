function [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, alpha0, n0, d, beta, theta)
% [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, nmax, alpha0, n0, d, beta, theta)
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
options=optimoptions(options, 'MaxFunEvals',10000);
options=optimoptions(options,'MaxIter',10000);
[x,fval,exitflag,output] = fmincon(@(x) TorqueNorm(kf, km, beta, theta, L, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconM(x, kf, km, beta, theta, L, d),  options);

% Substitution:
if exitflag == 1
    alphastar = [x(1) x(2) x(3) x(4)];
    nstar = [x(5) x(6) x(7) x(8)].';
else
    alphastar = alpha0;
    nstar = n0;
end


end
function [squareNormM] = TorqueNorm(kf, km, beta, theta, L, x)
M = Torque(kf, km, theta, beta, L, x);
% Objecive function maximize the squared norm of the torque M:
squareNormM = -(M(1)^2 + M(2)^2 + M(3)^2);
end
function [M] = Torque(kf, km, theta, beta, L, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
alpha = [x(1); x(2); x(3); x(4)];
n = [x(5); x(6); x(7); x(8)];
bRp1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
bRp2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
bRp3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
bRp4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
Op1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*[L 0 0].';
Op2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*[L 0 0].';
Op3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*[L 0 0].';
Op4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*[L 0 0].';
Tp1 = [0 0 kf*n(1)^2].'; % Thrust vector propeller 1
Tp2 = [0 0 kf*n(2)^2].'; % Thrust vector propeller 2
Tp3 = [0 0 kf*n(3)^2].'; % Thrust vector propeller 3
Tp4 = [0 0 kf*n(4)^2].'; % Thrust vector propeller 4
Tauext1 = [0 0 -km*n(1)^2].'; % Thrust vector propeller 1
Tauext2 = [0 0 km*n(2)^2].'; % Thrust vector propeller 2
Tauext3 = [0 0 -km*n(3)^2].'; % Thrust vector propeller 3
Tauext4 = [0 0 km*n(4)^2].'; % Thrust vector propeller 4
Taub = cross(Op1,bRp1*Tp1) + cross(Op2,bRp2*Tp2) + cross(Op3,bRp3*Tp3) + cross(Op4,bRp4*Tp4);
M = bRp1*Tauext1 + bRp2*Tauext2 + bRp3*Tauext3 + bRp4*Tauext4 + Taub;
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

