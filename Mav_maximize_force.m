function [alphastar, wstar] = Mav_maximize_force(n, kf, wmin, wmax, alphamin, alphamax, alpha0, w0, d, beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%MAV_MAXIMIZE_FORCE find optimal tilting angles and rotor speed that
%maximize force in direction d
%   Optimize alpha and n so the drone produce the maximal force in arbitrary direction d

%% Condition  Ax <= b        
A = []; 
b = [];                 

%% condition: Aeq.x = beq     
Aeq = [];
beq = [];

%% condition: lb <= x <= ub 
lb = zeros(2*n,1);
ub = zeros(2*n,1);
for i = n:-1:1
    lb(i) = alphamin; % lower bound
    lb(i+n) = wmin;
    ub(i) = alphamax; % upper bound
    ub(i+n) = wmax;
end

%% initial guess:
x0 = [alpha0; w0];

%% Objective function (maximize the 2-norm of the force produced by the propellers in the body frame):
% maximize norm of the force in an arbitrairy direction d:                                                                                                                                                               
% x = [alpha1, alpha2, alpha3, ..., n1, n2, n3, ....]
fun = @(x) -norm(thrust(x, n, kf, theta, beta, wRb));
            
%% optimization options
options = optimoptions('fmincon', 'Display', Display,'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

%% actual optimization
xstar = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, n, kf, theta, beta, d, wRb),  options);

%% Solution of the optimization
alphastar = xstar(1:n);
wstar = xstar(n+1:2*n);

%% Non linear constraints function 
function [c,ceq] = nonlinconF(x, n, kf, theta, beta, d, wRb)
% function [c,ceq] = nonlinconF(x, kf, theta, beta, d, wRb)
% defines the non linear constraints of the optimization problem

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Thrusts produced by the propellers in the body frame
T = thrust(x, n, kf, theta, beta, wRb);
  
% Constraint: c(x) <= 0
c = [];

% Constraint: ceq(x) = 0
% Thrust parallel to d conditions -> T x d = 0: 
ceq(1) = T(2)*d(3) - T(3)*d(2);
ceq(2) = T(3)*d(1) - T(1)*d(3);
ceq(3) = T(1)*d(2) - T(2)*d(1);
end

function [T] = thrust(x, n, kf, theta, beta, wRb)
% function [T] = thrust(x, kf, theta, beta, wRb)
% defines the thrust produced by the mav's propellers

% interval between arms in normal n-copter configuration
interval = 2*pi/n;
% pre-allocation:
bRp = zeros(3,3,n);
Tp = zeros(3,n);
T = [0; 0; 0];
alpha = zeros(1,n);
w = zeros(n,1);
for j =1:n
    alpha(j) = x(j);
    w(j) = x(n+j);
    %% Find the propellers rotation matrix
    bRp(:,:,j) = Rotz((j-1)*interval)*Rotz(theta(j))*Roty(beta(j))*Rotx(alpha(j));
    
    %% Find forces applied by all propellers thrusts on the body
    Tp(:,j) = [0 0 kf*w(j)^2].'; % force applied by every propeller in propeller frame
    T  = T + bRp(:,:,j)*Tp(:,j); % force applied by all the propellers in body frame
end
%% Thrust in the body frame
T = wRb*T;
end
end