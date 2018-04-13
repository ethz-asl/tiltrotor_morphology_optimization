function [alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter)
%[alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_MAX_THRUST find optimal tilting angles and Rotor
%speed

%% Optimization of alpha and n
% maximize norm of the Thrust F in an arbitrairy direction d:
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

options = optimoptions('fmincon', 'Display', Display,'Algorithm',Algorithm);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);
[x,fval,exitflag,output] = fmincon(@(x) minusThrustNorm(kf, theta, beta, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf, theta, beta, d),  options);
% verify that the found solution satisfies the constraints
if x(1)>=lb(1) && x(1)<=ub(1) && x(2)>=lb(2) && x(2)<=ub(2)  && x(3)>=lb(3) && x(3)<=ub(3) && x(4)>=lb(4) && x(4)<=ub(4) && x(5)>=lb(5) && x(5)<=ub(5) && x(6)>=lb(6) && x(6)<=ub(6)  && x(7)>=lb(7) && x(7)<=ub(7) && x(8)>=lb(8) && x(8)<=ub(8)
    alphastar = [x(1) x(2) x(3) x(4)];
    nstar = [x(5) x(6) x(7) x(8)].';
else
    alphastar = alpha0;
    nstar = n0;
end
end
function [minus_squareNormT] = minusThrustNorm(kf, beta, theta, x)
T = Thrust(kf, theta, beta, x);
% Objecive function maximize the squared norm of the thrust T:
minus_squareNormT = -(T(1)^2 + T(2)^2 + T(3)^2);
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
function [c,ceq] = nonlinconF(x, kf, theta, beta, d)
% function [c,ceq] = nonlincon(x, kf,d)
% % Hover condition
% if d(3) < 0
%     c = [];
% else
%     c(1) = -cos(x(1))*x(5)^2 - cos(x(2))*x(6)^2 - cos(x(3))*x(7)^2 - cos(x(4))*x(8)^2;
% end
% No hover cdt
c = [];
% Thrust parallel to d conditions: Ceq(x) = 0
T = Thrust(kf, theta, beta, x);
ceq(1) = T(2)*d(3) - T(3)*d(2);
ceq(2) = T(3)*d(1) - T(1)*d(3);
ceq(3) = T(1)*d(2) - T(2)*d(1);
end
