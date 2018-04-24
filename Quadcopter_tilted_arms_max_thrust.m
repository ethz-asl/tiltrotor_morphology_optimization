function [alphastar, nstar, exitflag] = Quadcopter_tilted_arms_max_thrust(kf, nmin, nmax, alphamin, alphamax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%[alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_MAX_THRUST find optimal tilting angles and Rotor speed
%   Optimize alpha and n so the drone produce the maximal force in arbitrary direction d

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
lb = [alphamin alphamin alphamin alphamin nmin nmin nmin nmin].'; % lower bound
ub = [alphamax alphamax alphamax alphamax nmax nmax nmax nmax].'; % upper bound

% initial guess:
x0 = [alpha0, n0.'];

% Objective function (maximize the 2-norm of the force produced by the propellers in the body frame):
fun = @(x) -sqrt((kf*(sin(x(1))*sin(theta(1)) + cos(x(1))*sin(beta(1))*cos(theta(1)))*x(5)^2 ...
                + kf*(sin(x(2))*cos(theta(2)) - cos(x(2))*sin(beta(2))*sin(theta(2)))*x(6)^2 ...
                - kf*(sin(x(3))*sin(theta(3)) + cos(x(3))*sin(beta(3))*cos(theta(3)))*x(7)^2 ...
                - kf*(sin(x(4))*cos(theta(4)) - cos(x(4))*sin(beta(4))*sin(theta(4)))*x(8)^2)^2 ...
           +(- kf*(sin(x(1))*cos(theta(1)) - cos(x(1))*sin(beta(1))*sin(theta(1)))*x(5)^2 ...
                + kf*(sin(x(2))*sin(theta(2)) + cos(x(2))*sin(beta(2))*cos(theta(2)))*x(6)^2 ...
                + kf*(sin(x(3))*(cos(theta(3))) - cos(x(3))*sin(beta(3))*sin(theta(3)))*x(7)^2 ...
                - kf*(sin(x(4))*sin(theta(4)) + cos(x(4))*sin(beta(4))*cos(theta(4)))*x(8)^2)^2 ...
            +(kf*cos(x(1))*cos(beta(1))*x(5)^2 + kf*cos(x(2))*cos(beta(2))*x(6)^2 ...
                + kf*cos(x(3))*cos(beta(3))*x(7)^2+ kf*cos(x(4))*cos(beta(4))*x(8)^2)^2);
            
% optimization options
options = optimoptions('fmincon', 'Display', Display,'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

% actual optimization
[xstar,fval,exitflag,output] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf, theta, beta, d),  options);

% Solution of the optimization
alphastar = [xstar(1) xstar(2) xstar(3) xstar(4)];
nstar = [xstar(5) xstar(6) xstar(7) xstar(8)].';

%% Tests
% T = Thrust(kf, theta, beta, xstar);
% lb = lb*1.01;
% ub = 1.01*ub;
% if xstar(1)>=lb(1) && xstar(1)<=ub(1) && xstar(2)>=lb(2) && xstar(2)<=ub(2)  && xstar(3)>=lb(3) && xstar(3)<=ub(3) && xstar(4)>=lb(4) && xstar(4)<=ub(4) && xstar(5)>=lb(5) && xstar(5)<=ub(5) && xstar(6)>=lb(6) && xstar(6)<=ub(6)  && xstar(7)>=lb(7) && xstar(7)<=ub(7) && xstar(8)>=lb(8) && xstar(8)<=ub(8)
% else
%     infoT = 'lb, ub not satisfied'
%     xstar = xstar
% end
% if ~isequal(round(cross(T,d)*10^2)/10^2,[0 0 0].')
%     infoT = '~//'
% end
% if norm([xstar(5), xstar(6), xstar(7), xstar(8)]) < 100
%     leiite = 1;
% end

%% Non linear constraints function 
function [c,ceq] = nonlinconF(x, kf, theta, beta, d)
% function [c,ceq] = nonlincon(x, kf,d)

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Thrusts produced by the propellers in the body frame
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
  
% Constraint: c(x) <= 0
c = [];

% Constraint: ceq(x) = 0
% Thrust parallel to d conditions -> T x d = 0: 
ceq(1) = T(2)*d(3) - T(3)*d(2);
ceq(2) = T(3)*d(1) - T(1)*d(3);
ceq(3) = T(1)*d(2) - T(2)*d(1);
end
end