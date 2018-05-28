function [nstar, betastar, thetastar, Lstar, obj_fun, exitflag] = Mav_optimize_n_beta_theta_L(dec, kf, km, nmin, nmax, Lmin, Lmax, g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, Display, ConstraintTolerance)
% [nstar, betastar, tetastar, Lstar, exitflag] = Mav_optimize_n_beta_theta_L(dec, kf, km, nmin, nmax, Lmin, Lmax, g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, Display, ConstraintTolerance)
%MAV_OPTIMIZE_BETA_THETA_L find optimal tilting angles and Rotor speed
%   Optimize alpha and n so the drone produce the maximal torque in arbitrary direction d

%% Optimization of alpha and n
% maximize norm of the Torque M in an arbitrairy direction d:

% x = [beta_1 beta_2 ... beta_n teta_1 teta_2 ... teta_4]

% Condition  Ax <= b
A = [];
b = [];

% Condition: Aeq.x = beq
Aeq = [];
beq = [];

%% condition: lb <= x <= ub
lb_n = nmin;
lb_beta = betamin*ones(1,nmax);
lb_theta = thetamin*ones(1,nmax);
lb_L = Lmin;
ub_n = nmax;
ub_beta = betamax*ones(1,nmax);
ub_theta = thetamax*ones(1,nmax);
ub_L = Lmax;

% fix the first arm position
lb_beta(1) = 0;
ub_beta(1) = 0;
lb_theta(1) = 0;
ub_theta(1) = 0;
% constraint the second arm to be on the same horizontal plan:
lb_beta(2) = 0;
ub_beta(2) = 0;

lb = [lb_n, lb_beta lb_theta lb_L];% lower bound
ub = [ub_n, ub_beta ub_theta ub_L];% upper bound

% optimization options
options = optimoptions('ga', 'Display', Display, 'ConstraintTolerance', ConstraintTolerance);

% actual optimization
[xstar, obj_fun, exitflag] = ga(@(x) objective_function(dec, x, kf, km, nmax, wmin, wmax, alphamin, alphamax, g, max_iterations), 2*nmax+2, A, b, Aeq, beq, lb, ub,[], 1,  options);
% [xstar,obj_fun,exitflag] = ga(@(x) objective_function(dec, n, x, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations),@(x) x(1), A,b,[],[], @(x) lb(x, nmin, betamin, thetamin, Lmin), @(x) ub(x, nmax, betamax, thetamax, Lmax),[], 1,options);
% Solution of the optimization
nstar = xstar(1);
betastar = xstar(2:nstar+1);
thetastar = xstar(nmax+2:nmax+1+nstar);
Lstar = xstar(2*nmax+2);

%% Objective function function
function [fun] = objective_function(dec, x, kf, km, nmax, wmin, wmax, alphamin, alphamax, g, max_iterations)
n = round(x(1));
beta = x(2:n+1);
theta = x(nmax+2:nmax+1+n);
L = x(2*nmax+2);

% Initial test to verify the consistence of the input:
size_beta = size(beta);
size_theta = size(theta);
if max(size_beta) ~= n &&  max(size_theta) ~= n
    fprintf('Arm angles defined not consistent with the number of arms')
    return;
end

% initialize the pitch yaw and roll angles to 0 (drone orientation w.r.t. to the world frame)
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
% Rotation Matrix mapping body frame to inertial frame
wRb = rotz(rad2deg(yaw0))*roty(rad2deg(pitch0))*rotz(rad2deg(roll0));

% Static matrix
% The static matrix are static allocation matrix that do not depend on the
% rotor orientation and speed.

% Vector containing all the decomposed vertical and horizontal forces:
% Fdec = [kf*cos(alpha(1))*w(1)^2; kf*sin(alpha(1))*w(1)^2;
%         kf*cos(alpha(2))*w(2)^2; kf*sin(alpha(2))*w(2)^2;
%         ...
%         kf*cos(alpha(n))*w(n)^2; kf*sin(alpha(n))*w(n)^2];

% The static matrix links Fdec to the force and the torque applied by the propellers to
% the drone body
% F = m*p'' = A_F_static*Fdec (w.r.t. to the body frame)
% M = Ib*wb' = A_M_static*Fdec (w.r.t. to the body frame)
[A_F_static, A_M_static] = Mav_static_matrix(kf, km, L, beta, theta, n, dec);

% The Moore-Penrose pseudo inverse of the static matrices allow to find
% Fdec from a desired force or torque applied on the drone.
% The rotor orientation and speed can then be deduced from Fdec
% => Fdec = inv(A_F_static)*Fdes
A_F_staticinv = pinv(A_F_static);
% => Fdec = inv(A_M_static)*Mdes
A_M_staticinv = pinv(A_M_static);

% Find the directions of the arms:
[m, ~, ~, ~, Op, ~] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
D =  [Op, -Op];
D_unit = D./vecnorm(D);
[~, length_D] = size(D_unit);

Fmin = []; % Matrix containing the maximum force appliable by the design in every direction of D
Mmin = []; % Matrix containing the maximum torques appliable by the design in every direction of D

% Find the maximum force and torque produced in the direction of the arms
for ii = 1:1:length_D
    d = D_unit(:,ii);
    % First, find the max force in direction d using static matrix
    Fdes = m*g*d;% Set the initial desired force to be the force to hover in direction d
    k=4; % Start with big steps
    for i = 1:max_iterations % Loop to find the maximal force appliable by the drone in direction d
        
        Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
        
        % Retrieve rotors speeds and orientations from Fdec
        [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
        
        alpha_bound = [alpha0(alpha0>alphamax), alpha0(alpha0<alphamin)]; % fill alpha_bound if the bounds of alpha are violated
        w_bound = [w0(w0>wmax), w0(w0<wmin)];% fill w_bound if the bounds of w are violated
        
        if isempty(w_bound) && isempty(alpha_bound) % if constraint respected
            % Slowly increase Fdes until the obtained alpha0 and w0 does
            % not respect their bounds anymore.
            Fdes = Fdes + k*d*(n*wmax^2*kf-m*g)/max_iterations;
        else% If alpha0 and w0 does not respect their bounds anymore.
            
            % Return to the previous Fdes
            Fdes = Fdes - k*d*(n*wmax^2*kf-m*g)/max_iterations;
            
            if k < 0.25 % if step size under a treshold break.
                break;
            else
                k = k/2; % else diminish the step size and loop agai.
            end
        end
    end
    Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
    
    % Retrieve rotors speeds and orientations from Fdec
    [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
    
    %calculate linear acceleration with this alphastar and nstar
    [~, ~,pdotdot, ~] = Mav_dynamic(n, kf, km, wRb, alpha0, beta, theta,w0, L, g, dec, false);
    F0 = m*pdotdot;
    Fmin(ii) = norm(F0);
    
    % find max torque in direction d using static matrix
    % Set the initial desired torque to be the torque produced if the force
    % to hover was applied at the end of one of the arms
    Mdes = d*(m*g*L);
    % Loop to find the maximal torque appliable by the drone in direction d
    k=4;
    for i = 1:max_iterations
        Fdec = A_M_staticinv*Mdes; % Fdec = inv(Astatic)*Fdes
        
        % Retrieve rotors speeds and orientations from Fdec
        [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
        
        alpha_bound = [alpha0(alpha0>alphamax), alpha0(alpha0<alphamin)];
        w_bound = [w0(w0>wmax), w0(w0<wmin)];
        if isempty(w_bound) && isempty(alpha_bound)
            % Slowly increase Fdes until the obtained alpha0 and w0 does
            % not respect their bounds anymore.
            Mdes = Mdes + k*d*(n*L*wmax^2*kf-m*g*L)/max_iterations;
        else
            % If alpha0 and w0 does not respect their bounds anymore.
            % Return to the previous Fdes
            Mdes = Mdes - k*d*(n*L*wmax^2*kf-m*g*L)/max_iterations;
            if k < 0.25
                break;
            else
                k = k/2;
            end
        end
    end
    Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
    
    % Retrieve rotors speeds and orientations from Fdec
    [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
    
    % calculate angular acceleration with this alphastar and nstar
    [~, Ib, ~, wbdot] = Mav_dynamic(n, kf, km, wRb, alpha0, beta, theta, w0, L, g, dec, false);
    wbdot = round(dec*wbdot)/dec;
    M0 = Ib*wbdot;
    Mmin(ii) = norm(M0);
end
% Compute the inertia as a fct of L, beta and theta
[m, Ib] = Mav_inertias(n, L, theta, beta);
%% Objective function
fun = -sum(Mmin) -sum(Fmin) + 600*sqrt(Ib(1,1)^2 + Ib(2,2)^2 + Ib(3,3)^2);
end
end