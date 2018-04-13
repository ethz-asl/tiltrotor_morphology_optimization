function [alphastar, nstar] = Quadcopter_tilted_arms_min_hover(kf, m, g, fdes, nmax, alpha0, n0, beta, theta, Display, Algorithm, maxIter)
% [alphastar, nstar] = Quadcopter_tilted_arms_opt_hover(kf,m, g, nmax, alpha0, n0, beta, theta, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_MIN_HOVER finds the optimal hover mode in direction d (parallel to fdes)
%   Optimize alpha and n to obtain the most efficient hover when drone is
%   oriented in direction z parallel to fdes.
%% Optimization of alpha and n
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
% optimization options 
options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm);
options=optimoptions(options, 'MaxFunEvals',maxIter);
options=optimoptions(options,'MaxIter',maxIter);
% actual optimization 
xstar = fmincon(@(x) minRotorSpeed(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlincon(x, kf, theta, beta, fdes), options);

% verify that the found solution satisfies the constraints
T = Thrust(kf, theta, beta, xstar);
if xstar(1)>=lb(1) && xstar(1)<=ub(1) && xstar(2)>=lb(2) && xstar(2)<=ub(2)  && xstar(3)>=lb(3) && xstar(3)<=ub(3) && xstar(4)>=lb(4) && xstar(4)<=ub(4) && xstar(5)>=lb(5) && xstar(5)<=ub(5) && xstar(6)>=lb(6) && xstar(6)<=ub(6)  && xstar(7)>=lb(7) && xstar(7)<=ub(7) && xstar(8)>=lb(8) && xstar(8)<=ub(8) && norm(T)>= norm(fdes)  && isequal(round(cross(T,fdes)*10^2)/10^2,[0 0 0].')
    alphastar = [xstar(1) xstar(2) xstar(3) xstar(4)];
    nstar = [xstar(5) xstar(6) xstar(7) xstar(8)].';
else
    alphastar = alpha0;
    nstar = n0;
end

end
function [minN] = minRotorSpeed(x)
% Objecive function minimize the rotor speeds 
minN = x(5) + x(6) + x(7) + x(8);
end
function [T] = Thrust(kf, theta, beta, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
alpha = [x(1); x(2); x(3); x(4)];
n = [x(5); x(6); x(7); x(8)];
%% Thrusts of the propellers in the body frame
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
% function [c,ceq] = nonlincon(x, kf,d)
T = Thrust(kf, theta, beta, x);
% Condition  c(x) <= 0
c = [];
% Condition  ceq(x) = 0
% force applied by the propellers equal to the force needed by the drone to hover
ceq(1) = T(1)-fdes(1);
ceq(2) = T(2)-fdes(2);
ceq(3) = T(3)-fdes(3);
end


