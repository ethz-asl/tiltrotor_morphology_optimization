function [alphastar, nstar, exitflag, T] = Quadcopter_tilted_arms_min_hover(kf, fdes, nmin, nmax, alphamin, alphamax, alpha0, n0, beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, min_alpha)
% [alphastar, nstar, exitflag, T] = Quadcopter_tilted_arms_min_hover(kf, m, g, fdes, nmin, nmax, alphamin, alphamax, alpha0, n0, beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, min_alpha)
%QUADCOPTER_TILTED_ARMS_MIN_HOVER finds the optimal hover mode in direction d (parallel to fdes)
%   Optimize alpha and n to obtain the most efficient hover when drone is
%   oriented in direction z parallel to fdes.

%% Optimization of alpha and n

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Condition  Ax <= b
A = [];
b = [];

% Condition: Aeq.x = beq
Aeq = [];
beq = [];

% Condition: lb <= x <= ub 
lb = [alphamin alphamin alphamin alphamin nmin nmin nmin nmin].'; % lower bound
ub = [alphamax alphamax alphamax alphamax nmax nmax nmax nmax].'; % upper bound

% Initial guess:
x0 = [alpha0, n0.'];

% Objective function:
if min_alpha % only minimise the tilting angle
    fun = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2); % minimize 2-norm of alpha
else % only minimise the proppeler speed
    fun = @(x) sqrt(x(5)^2 + x(6)^2 + x(7)^2 + x(8)^2); % minimize 2-norm of n
end

% optimization options
options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals',maxIter);
options=optimoptions(options,'MaxIter',maxIter);

% actual optimization
[x,fval,exitflag,output] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlincon(x, kf, theta, beta, fdes), options);

% Thrusts of the propellers in the body frame after optimisation
T = [(kf*(sin(x(1))*sin(theta(1)) + cos(x(1))*sin(beta(1))*cos(theta(1)))*x(5)^2 ...
      +kf*(sin(x(2))*cos(theta(2)) - cos(x(2))*sin(beta(2))*sin(theta(2)))*x(6)^2 ...
      -kf*(sin(x(3))*sin(theta(3)) + cos(x(3))*sin(beta(3))*cos(theta(3)))*x(7)^2 ...
      -kf*(sin(x(4))*cos(theta(4)) - cos(x(4))*sin(beta(4))*sin(theta(4)))*x(8)^2); ...
     (-kf*(sin(x(1))*cos(theta(1)) - cos(x(1))*sin(beta(1))*sin(theta(1)))*x(5)^2 ...
      +kf*(sin(x(2))*sin(theta(2)) + cos(x(2))*sin(beta(2))*cos(theta(2)))*x(6)^2 ...
      +kf*(sin(x(3))*(cos(theta(3))) - cos(x(3))*sin(beta(3))*sin(theta(3)))*x(7)^2 ...
      -kf*(sin(x(4))*sin(theta(4)) + cos(x(4))*sin(beta(4))*cos(theta(4)))*x(8)^2); ...
     (kf*cos(x(1))*cos(beta(1))*x(5)^2 + kf*cos(x(2))*cos(beta(2))*x(6)^2 ...
      +kf*cos(x(3))*cos(beta(3))*x(7)^2+ kf*cos(x(4))*cos(beta(4))*x(8)^2)];

% Solution of the optimization
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';

%% Tests
% if x(1)>=lb(1) && x(1)<=ub(1) && x(2)>=lb(2) && x(2)<=ub(2)  && x(3)>=lb(3) && x(3)<=ub(3) && x(4)>=lb(4) && x(4)<=ub(4) && x(5)>=lb(5) && x(5)<=ub(5) && x(6)>=lb(6) && x(6)<=ub(6)  && x(7)>=lb(7) && x(7)<=ub(7) && x(8)>=lb(8) && x(8)<=ub(8)
% else
%     infoH = 'lb, ub not satisfied'
%     xstar = x
% end
% if ~isequal(round(T*10^2)/10^2,round(fdes*10^2)/10^2)
%     infoH = '~//'
% end

%% Non linear constraints function 
function [c,ceq] = nonlincon(x, kf, theta, beta, fdes)
% function [c,ceq] = nonlincon(x, kf, theta, beta, fdes, n0)

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Thrusts of the propellers in the body frame
T = [(kf*(sin(x(1))*sin(theta(1)) + cos(x(1))*sin(beta(1))*cos(theta(1)))*x(5)^2 ...
      +kf*(sin(x(2))*cos(theta(2)) - cos(x(2))*sin(beta(2))*sin(theta(2)))*x(6)^2 ...
      -kf*(sin(x(3))*sin(theta(3)) + cos(x(3))*sin(beta(3))*cos(theta(3)))*x(7)^2 ...
      -kf*(sin(x(4))*cos(theta(4)) - cos(x(4))*sin(beta(4))*sin(theta(4)))*x(8)^2); ...
     (-kf*(sin(x(1))*cos(theta(1)) - cos(x(1))*sin(beta(1))*sin(theta(1)))*x(5)^2 ...
      +kf*(sin(x(2))*sin(theta(2)) + cos(x(2))*sin(beta(2))*cos(theta(2)))*x(6)^2 ...
      +kf*(sin(x(3))*(cos(theta(3))) - cos(x(3))*sin(beta(3))*sin(theta(3)))*x(7)^2 ...
      -kf*(sin(x(4))*sin(theta(4)) + cos(x(4))*sin(beta(4))*cos(theta(4)))*x(8)^2); ...
     (kf*cos(x(1))*cos(beta(1))*x(5)^2 + kf*cos(x(2))*cos(beta(2))*x(6)^2 ...
      +kf*cos(x(3))*cos(beta(3))*x(7)^2+ kf*cos(x(4))*cos(beta(4))*x(8)^2)];

% Condition  c(x) <= 0
c = [];

% Condition:  ceq(x) = 0

% Thrust parallel to fdes conditions -> T x fdes = 0: 
ceq(1) = T(2)*fdes(3) - T(3)*fdes(2);
ceq(2) = T(3)*fdes(1) - T(1)*fdes(3);
ceq(3) = T(1)*fdes(2) - T(2)*fdes(1);

% Thrust equal to fdes conditions -> fdes - T = 0: 
ceq(4) = fdes(1)-T(1);
ceq(5) = fdes(2)-T(2);
ceq(6) = fdes(3)-T(3);
% ceq = [];
end
end
