function [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter)
% [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter)
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

options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);
[x,fval,exitflag,output] = fmincon(@(x) TorqueNorm(kf, km, beta, theta, L, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconM(x, kf, km, beta, theta, L, d),  options);
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
% % Substitution:
% if exitflag == 1
%     alphastar = [x(1) x(2) x(3) x(4)];
%     nstar = [x(5) x(6) x(7) x(8)].';
% else
%     alphastar = alpha0;
%     nstar = n0;
% end


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
%% Time consuming methode:
% bRp1 = Rotz(rad2deg(theta(1)))*Roty(rad2deg(beta(1)))*Rotx(rad2deg(alpha(1)));
% bRp2 = Rotz(rad2deg(pi/2+theta(2)))*Roty(rad2deg(beta(2)))*Rotx(rad2deg(alpha(2)));
% bRp3 = Rotz(rad2deg(pi+theta(3)))*Roty(rad2deg(beta(3)))*Rotx(rad2deg(alpha(3)));
% bRp4 = Rotz(rad2deg(3*pi/2+theta(4)))*Roty(rad2deg(beta(4)))*Rotx(rad2deg(alpha(4)));
% Op1 = Rotz(rad2deg(theta(1)))*Roty(rad2deg(beta(1)))*[L 0 0].';
% Op2 = Rotz(rad2deg(pi/2+theta(2)))*Roty(rad2deg(beta(2)))*[L 0 0].';
% Op3 = Rotz(rad2deg(pi+theta(3)))*Roty(rad2deg(beta(3)))*[L 0 0].';
% Op4 = Rotz(rad2deg(3*pi/2+theta(4)))*Roty(rad2deg(beta(4)))*[L 0 0].';
% Tp1 = [0 0 kf*n(1)^2].'; % Thrust vector propeller 1
% Tp2 = [0 0 kf*n(2)^2].'; % Thrust vector propeller 2
% Tp3 = [0 0 kf*n(3)^2].'; % Thrust vector propeller 3
% Tp4 = [0 0 kf*n(4)^2].'; % Thrust vector propeller 4
% Tauext1 = [0 0 -km*n(1)^2].'; % Thrust vector propeller 1
% Tauext2 = [0 0 km*n(2)^2].'; % Thrust vector propeller 2
% Tauext3 = [0 0 -km*n(3)^2].'; % Thrust vector propeller 3
% Tauext4 = [0 0 km*n(4)^2].'; % Thrust vector propeller 4
% Taub = cross(Op1,bRp1*Tp1) + cross(Op2,bRp2*Tp2) + cross(Op3,bRp3*Tp3) + cross(Op4,bRp4*Tp4);
% M = bRp1*Tauext1 + bRp2*Tauext2 + bRp3*Tauext3 + bRp4*Tauext4 + Taub;
%% Pre-computed method (Faster)
Taub1 = - L*kf*n(1)^2*(sin(beta(1))*sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(theta(1))) ...
          + L*kf*n(2)^2*(sin(beta(2))*sin(alpha(2))* sin(theta(2)) + cos(alpha(2))*cos(theta(2))) ...
          + L*kf*n(3)^2*(sin(beta(3))*sin(alpha(3))*cos(theta(3)) - cos(alpha(3))* sin(theta(3))) ...
          - L*kf*n(4)^2*(sin(beta(4))*sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*cos(theta(4)));
Taub2 = - L*kf*n(1)^2*(sin(beta(1))*sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*cos(theta(1))) ...
          - L*kf*n(2)^2*(sin(beta(2))*sin(alpha(2))*cos(theta(2)) - cos(alpha(2))* sin(theta(2))) ...
          + L*kf*n(3)^2*(sin(beta(3))*sin(alpha(3))* sin(theta(3)) + cos(alpha(3))*cos(theta(3))) ...
          + L*kf*n(4)^2*(sin(beta(4))*sin(alpha(4))*cos(theta(4)) - cos(alpha(4))*sin(theta(4)));
Taub3 = - L*kf*n(1)^2*cos(beta(1))*sin(alpha(1)) ...
          - L*kf*n(2)^2*cos(beta(2))*sin(alpha(2)) ...
          - L*kf*n(3)^2*cos(beta(3))*sin(alpha(3)) ...
          - L*kf*n(4)^2*cos(beta(4))*sin(alpha(4));
M1 = -km*(sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*sin(beta(1))*cos(theta(1)))*n(1)^2 ...
       + km*(sin(alpha(2))*(cos(theta(2))) - cos(alpha(2))*sin(beta(2))*sin(theta(2)))*n(2)^2 ...
       + km*(sin(alpha(3))*sin(theta(3)) + cos(alpha(3))*sin(beta(3))*cos(theta(3)))*n(3)^2 ...
       - km*(sin(alpha(4))*(cos(theta(4))) - cos(alpha(4))*sin(beta(4))*sin(theta(4)))*n(4)^2;
M2 = + km*(sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(beta(1))*sin(theta(1)))*n(1)^2 ...
     + km*(sin(alpha(2))*sin(theta(2)) + cos(alpha(2))*sin(beta(2))*(cos(theta(2))))*n(2)^2 ...
       - km*(sin(alpha(3))*cos(theta(3)) - cos(alpha(3))*sin(beta(3))*sin(theta(3)))*n(3)^2 ...
       - km*(sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*sin(beta(4))*cos(theta(4)))*n(4)^2;
M3 = -km*cos(alpha(1))*cos(beta(1))*n(1)^2 + km*cos(alpha(2))*cos(beta(2))*n(2)^2 ...
       - km*cos(alpha(3))*cos(beta(3))*n(3)^2 + km*cos(alpha(4))*cos(beta(4))*n(4)^2;
M = [M1+Taub1; M2+Taub2; M3+Taub3];
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

