function [alphastar, nstar,i] = Quadcopter_tilted_arms_min_hover(i,kf, m, g, fdes, nmax, alpha0, n0, beta, theta, Display, Algorithm, maxIter)
% [alphastar, nstar] = Quadcopter_tilted_arms_opt_hover(kf,m, g, nmax, alpha0, n0, beta, theta, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_MIN_HOVER finds the optimal hover mode 
%   Optimize alpha and n to obtain the most efficient hover when drone is
%   oriented in direction d
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
options=optimoptions(options, 'MaxFunEvals',maxIter);
options=optimoptions(options,'MaxIter',maxIter);
[x,fval,exitflag,output] = fmincon(@(x) ThrustNorm(kf, fdes, theta, beta, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlincon(x, kf, theta, beta, fdes), options);
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
% Substitution:
% T = Thrust(kf, theta, beta, x);
% if norm(T) < 0.99*m*g
%     i = i+1;
%     T
%     fdes
% %     alphastar = alpha0;
% %     nstar = n0;
% end

if exitflag == 1
    i = i+1;
    alphastar = [x(1) x(2) x(3) x(4)];
    nstar = [x(5) x(6) x(7) x(8)].';
else
    alphastar = alpha0;
    nstar = n0;
end


end
function [squareNormT] = ThrustNorm(kf, fdes, beta, theta, x)
T = Thrust(kf, theta, beta, x);
% Objecive function maximize the squared norm of the thrust T:
squareNormT = (T(1))^2 + (T(2))^2 + (T(3))^2;

end
function [T] = Thrust(kf, theta, beta, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
alpha = [x(1); x(2); x(3); x(4)];
n = [x(5); x(6); x(7); x(8)];
%% Time consuming methode:
% bRp1 = Rotz(rad2deg(theta(1)))*Roty(rad2deg(beta(1)))*Rotx(rad2deg(alpha(1)));
% bRp2 = Rotz(rad2deg(pi/2+theta(2)))*Roty(rad2deg(beta(2)))*Rotx(rad2deg(alpha(2)));
% bRp3 = Rotz(rad2deg(pi+theta(3)))*Roty(rad2deg(beta(3)))*Rotx(rad2deg(alpha(3)));
% bRp4 = Rotz(rad2deg(3*pi/2+theta(4)))*Roty(rad2deg(beta(4)))*Rotx(rad2deg(alpha(4)));
% Tp1 = [0 0 kf*n(1)^2].'; % Thrust vector propeller 1
% Tp2 = [0 0 kf*n(2)^2].'; % Thrust vector propeller 2
% Tp3 = [0 0 kf*n(3)^2].'; % Thrust vector propeller 3
% Tp4 = [0 0 kf*n(4)^2].'; % Thrust vector propeller 4
% T = (bRp1*Tp1 + bRp2*Tp2 + bRp3*Tp3 + bRp4*Tp4);
%% Pre-computed method (Faster)
T1 = kf*(sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*sin(beta(1))*cos(theta(1)))*n(1)^2 ...
       + kf*(sin(alpha(2))*cos(theta(2)) - cos(alpha(2))*sin(beta(2))*sin(theta(2)))*n(2)^2 ...
       - kf*(sin(alpha(3))*sin(theta(3)) + cos(alpha(3))*sin(beta(3))*cos(theta(3)))*n(3)^2 ...
       - kf*(sin(alpha(4))*cos(theta(4)) - cos(alpha(4))*sin(beta(4))*sin(theta(4)))*n(4)^2;

T2 = - kf*(sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(beta(1))*sin(theta(1)))*n(1)^2 ...
       + kf*(sin(alpha(2))*sin(theta(2)) + cos(alpha(2))*sin(beta(2))*cos(theta(2)))*n(2)^2 ...
       + kf*(sin(alpha(3))*(cos(theta(3))) - cos(alpha(3))*sin(beta(3))*sin(theta(3)))*n(3)^2 ...
       - kf*(sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*sin(beta(4))*cos(theta(4)))*n(4)^2;
T3 = kf*cos(alpha(1))*cos(beta(1))*n(1)^2 + kf*cos(alpha(2))*cos(beta(2))*n(2)^2 ...
       + kf*cos(alpha(3))*cos(beta(3))*n(3)^2+ kf*cos(alpha(4))*cos(beta(4))*n(4)^2;
T = [T1;T2;T3];
end
function [c,ceq] = nonlincon(x, kf, theta, beta, fdes)
T = Thrust(kf, theta, beta, x);
% function [c,ceq] = nonlincon(x, kf,d)
c = [];
c(1) = -T(1) - T(2) -T(3)+fdes(1) + fdes(2) +fdes(3);
% Thrust equal to fdes conditions: Ceq(x) = 0
ceq = [];
ceq(1) = T(1)-fdes(1);
ceq(2) = T(2)-fdes(2);
ceq(3) = T(3)-fdes(3);
% ceq(4) = T(2)*fdes(3) - T(3)*fdes(2);
% ceq(5) = T(3)*fdes(1) - T(1)*fdes(3);
% ceq(6) = T(1)*fdes(2) - T(2)*fdes(1);

end


