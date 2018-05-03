function [alphastar,wstar] = Mav_maximize_torque(n, kf, km, L, wmin, wmax, alphamin, alphamax, alpha0, w0, d, beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%MAV_MAXIMIZE_TORQUE find optimal tilting angles and rotor speed that
%maximize torque in direction d
%   Optimize alpha and n so the drone produce the maximal torque in arbitrary direction d

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

%% Objective function (maximize the 2-norm of the torque produced by the propellers in the body frame):
% maximize norm of the torque in an arbitrairy direction d:                                                                                                                                                               
% x = [alpha1, alpha2, alpha3, ..., n1, n2, n3, ....]
fun = @(x) -norm(torque(x, n, L, km, kf, theta, beta));
            
%% optimization options
options = optimoptions('fmincon', 'Display', Display,'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

%% actual optimization
xstar = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, n, L, km, kf, theta, beta, d),  options);

%% Solution of the optimization
alphastar = xstar(1:n);
wstar = xstar(n+1:2*n);

%% Non linear constraints function 
function [c,ceq] = nonlinconF(x, n, L, km, kf, theta, beta, d)
% function [c,ceq] = nonlinconF(x, kf, theta, beta, d)
% defines the non linear constraints of the optimization problem

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Thrusts produced by the propellers in the body frame
M = torque(x, n, L, km, kf, theta, beta);
  
% Constraint: c(x) <= 0
c = [];

% Constraint: ceq(x) = 0
% Torque parallel to d conditions -> M x d = 0: 
ceq(1) = M(2)*d(3) - M(3)*d(2);
ceq(2) = M(3)*d(1) - M(1)*d(3);
ceq(3) = M(1)*d(2) - M(2)*d(1);
end

function [M] = torque(x, n, L, km, kf, theta, beta)
% function [M] = torque(x, kf, theta, beta, wRb)
% defines the torque produced by the mav's propellers

% interval between arms in normal n-copter configuration
interval = 2*pi/n;
% pre-allocation:
bRp = zeros(3,3,n);
Op = zeros(3,n);
Tp = zeros(3,n);
Tauext = zeros(3,n);
M = [0; 0; 0];
alpha = zeros(1,n);
w = zeros(n,1);
for j =1:n
    alpha(j) = x(j);
    w(j) = x(n+j);
    %% Find the propellers rotation matrix
    bRp(:,:,j) = Rotz((j-1)*interval)*Rotz(theta(j))*Roty(beta(j))*Rotx(alpha(j));
    %% Find the propellers positions in the body frame
    Op(:,j) = Rotz((j-1)*interval)*Rotz(theta(j))*Roty(beta(j))*[L 0 0].';
    %% Find forces applied by all propellers thrusts on the body
    Tp(:,j) = [0 0 kf*w(j)^2].'; % force applied by every propeller in propeller frame
    %% Find torques applied by all propellers on the body (in propeller frame)
    if mod(j,2) == 0
        c = 1; % counter torque defined negative for clock wise rotating propellers (all odd propellers)
    else
        c = -1;
    end
    Tauext(:,j) = [0 0 c*km*w(j)^2].'; % counter torque produced by every propeller in propeller frame
    %% Torque in the body frame
    M = M + bRp(:,:,j)*Tauext(:,j) + cross(Op(:,j),bRp(:,:,j)*Tp(:,j)); % torque applied by every propeller on the body
end
end
end
