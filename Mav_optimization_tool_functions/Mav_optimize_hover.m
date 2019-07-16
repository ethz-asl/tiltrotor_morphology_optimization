function [alphastar, wstar, exitflag, TH] = Mav_optimize_hover(n, kf, fdes, wmin, wmax, alphamin, alphamax, alpha0, w0, beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%MAV_OPTIMIZE_HOVER find optimal tilting angles and rotor speed that
%maximize the hover efficiency in direction d

%% Condition  Ax <= b        
A = []; 
b = [];                 

%% condition: Aeq.x = beq     
Aeq = [];
beq = [];

%% condition: lb <= x <= ub 
lb_alpha = alphamin*ones(n,1);
ub_alpha = alphamax*ones(n,1);
lb_w = wmin*ones(n,1);
ub_w = wmax*ones(n,1);
lb = [lb_alpha; lb_w];
ub = [ub_alpha; ub_w];

%% initial guess:
x0 = [alpha0; w0];

%% Objective function (minimize the 2-norm of the propeller rotation speed):
% minimize the norm of the propeller speeds in order to have an energy efficient hover                                                                                                                                                             
% x = [alpha1, alpha2, alpha3, ..., n1, n2, n3, ....]

% fun = @(x) norm(x(n+1:2*n));

%% optimization options
options = optimoptions('fmincon', 'Display', Display,'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

%% actual optimization
[xstar,~,exitflag,~]  = fmincon(@(x) objective_function(n,x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, n, kf, theta, beta, fdes, wRb),  options);

%% Solution of the optimization
alphastar = xstar(1:n);
wstar = xstar(n+1:2*n);
TH = thrust(xstar, n, kf, theta, beta, wRb);

%% Non linear constraints function 
function [c,ceq] = nonlinconF(x, n, kf, theta, beta, fdes, wRb)
% function [c,ceq] = nonlinconF(x, kf, theta, beta, d, wRb)
% defines the non linear constraints of the optimization problem

% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]

% Thrusts produced by the propellers in the body frame
T = thrust(x, n, kf, theta, beta, wRb);
  
% Constraint: c(x) <= 0
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
end

function [fun] = objective_function(n,x)
%% Objective function (minimize the 2-norm of the propeller rotation speed):
% minimize the norm of the propeller speeds in order to have an energy efficient hover                                                                                                                                                             
% x = [alpha1, alpha2, alpha3, ..., n1, n2, n3, ....]
fun = 0;
for i = n+1:2*n
    fun = fun + x(i)^2;
end
fun = sqrt(fun);
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
alpha = x(1:n);
w = x(n+1:2*n);
for l =1:n
    %% Find the propellers rotation matrix
    bRp(:,:,l) = Rotz((l-1)*interval)*Rotz(theta(l))*Roty(beta(l))*Rotx(alpha(l));
    
    %% Find forces applied by all propellers thrusts on the body
    Tp(:,l) = [0 0 kf*w(l)^2].'; % force applied by every propeller in propeller frame
    T  = T + bRp(:,:,l)*Tp(:,l); % force applied by all the propellers in body frame
end
%% Thrust in the body frame
T = wRb*T;
end
end

