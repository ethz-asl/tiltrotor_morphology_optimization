function [betastar, thetastar, Lstar, obj_fun, exitflag] = Mav_optimize_beta_theta_L(cost_fct_case, Optimize_theta, Optimize_L, direction, dec, n, kf, km, Lmin, Lmax, L0, g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
% [betastar, tetastar, Lstar, exitflag] = Mav_optimize_beta_theta_L(dec, n, kf, km, L, g, nmin, nmax, betamin, betamax, thetamin, thetamax, alphadotmax, alphamin, alphamax, opt_iterations, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
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
lb_beta = betamin*ones(1,n);
ub_beta = betamax*ones(1,n);
lb_theta = thetamin*ones(1,n);
ub_theta = thetamax*ones(1,n);
lb_L = Lmin;
ub_L = Lmax;
lb = lb_beta;% lower bound
ub = ub_beta;% upper bound
x0 = beta0; % initial guess
if Optimize_theta    
    lb = [lb lb_theta];% lower bound
    ub = [ub ub_theta];% upper bound
    x0 = [x0 theta0]; % initial guess
end
if Optimize_L   
    lb = [lb lb_L];% lower bound
    ub = [ub ub_L];% upper bound
    x0 = [x0 L0]; % initial guess
end
%% optimization options
options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

%% actual optimization
if cost_fct_case == '2'
     % Maximize Mmin and Fmin and minimize the inertia
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_Mmin_Fmin_min_I(Optimize_theta, Optimize_L, dec, n, x, L0, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations), x0, A, b, Aeq, beq, lb, ub,[],  options);
elseif cost_fct_case == '3'
    % Maximize F, M and Heff in every direction
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_every_direction(Optimize_theta, Optimize_L, dec, n, x, L0, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance), x0, A, b, Aeq, beq, lb, ub,[],  options);
elseif cost_fct_case == '4'
    % Maximize Fmax and Mmax in d direction
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_Mmax_Fmax(Optimize_theta, Optimize_L, dec, n, x, L0, direction, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations), x0, A, b, Aeq, beq, lb, ub,[],  options);
elseif cost_fct_case == '5'
    % Maximize Fmax and Mmax in d direction while hovering in every directions
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_Mmax_Fmax(Optimize_theta, Optimize_L, dec, n, x, L0, direction, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(Optimize_theta,Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance),  options);
elseif cost_fct_case == '6'
    % Maximize the force and torque in x, y and z directions
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_M_F_in_xyz(Optimize_theta, Optimize_L, dec, n, x, L0, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations), x0, A, b, Aeq, beq, lb, ub,[],  options);
elseif cost_fct_case == '7'
    % Maximize the force and torque space volumes.
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_F_M_vol(Optimize_theta, Optimize_L, dec, n, x, L0, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance), x0, A, b, Aeq, beq, lb, ub,[],  options);
else
    % Maximize Mmin and Fmin
    [xstar, obj_fun, exitflag, ~] = fmincon(@ (x) objective_function_max_Mmin_Fmin(Optimize_theta, Optimize_L, dec, n, x, L0, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations), x0, A, b, Aeq, beq, lb, ub,[],  options);
end

%% Solution of the optimization

Lstar = L0;
thetastar = theta0;
betastar = xstar(1:n);
if Optimize_theta 
    thetastar = xstar(n+1:2*n);
    if Optimize_L 
        Lstar = xstar(2*n+1);
    end
else
    if Optimize_L
        Lstar  = xstar(n+1);
    end
end
%% Objective function functions
%% Maximize Mmin and Fmin and minimize the inertia
function [fun] = objective_function_max_Mmin_Fmin(Optimize_theta, Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

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
        	if i == 1
               % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Fdes = 0*d;
               break;
            end
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
        else % If alpha0 and w0 does not respect their bounds anymore.
            if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Mdes = 0*d;
               break;
            end
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

%% Objecticve function
fun = -sum(Mmin) -sum(Fmin);
end

%% Maximize Mmin and Fmin and minimize the inertia
function [fun] = objective_function_max_Mmin_Fmin_min_I(Optimize_theta, Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

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
        	if i == 1% if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Fdes = 0*d;
               break;
            end
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
        else% If alpha0 and w0 does not respect their bounds anymore.
            if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Mdes = 0*d;
               break;
            end
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
[~, Ib] = Mav_inertias(n, L, theta, beta);
% fun1 = 2*sum(Mmin)
% fun2 = sum(Fmin)
%% Objecticve function
fun = -sum(Mmin) -sum(Fmin) + 300*norm(vecnorm(Ib));
end

%% Maximize Fmax and Mmax in d direction and minimize the inertia
function [fun] = objective_function_max_Mmax_Fmax(Optimize_theta, Optimize_L, dec, n, x, L, d, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

d = d./norm(d); % normalize d

%% Initial test to verify the consistence of the input:
size_beta = size(beta);
size_theta = size(theta);
if max(size_beta) ~= n &&  max(size_theta) ~= n
    fprintf('Arm angles defined not consistent with the number of arms')
    return;
end

%% initialize the pitch yaw and roll angles to 0 (drone orientation w.r.t. to the world frame)
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
% Rotation Matrix mapping body frame to inertial frame
wRb = rotz(rad2deg(yaw0))*roty(rad2deg(pitch0))*rotz(rad2deg(roll0));

%% Static matrix
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


%% Find the maximum force and torque produced in direction d 
[m, ~] = Mav_inertias(n, L, theta, beta);
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
     	if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Fdes = 0*d;
               break;
            end
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
Fz = norm(F0);

% find max torque in z direction using static matrix
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
    else % If alpha0 and w0 does not respect their bounds anymore.
      	if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Mdes = 0*d;
               break;
            end
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
Mz = norm(M0);

%% Objective function:
fun = -Mz -Fz;
end
%% Maximize Fmin, Mmin and Hover efficiency in every direction 
function [fun] = objective_function_every_direction(Optimize_theta, Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

%% Initial test to verify the consistence of the input:
size_beta = size(beta);
size_theta = size(theta);
if max(size_beta) ~= n &&  max(size_theta) ~= n
    fprintf('Arm angles defined not consistent with the number of arms')
    return;
end

%% Compute different metrix for a drone design
step = 0.25;
[~, ~, Heff, ~, ~, F,~, ~, ~, M] = Mav_compute_metrics(dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%% Objective function:
% fun = -F_vol -M_vol;
fun = -sum(vecnorm(F)) -sum(vecnorm(M)) -sum(Heff);
end
%% Maximize the volume of the force space and the torque space.
function [fun] = objective_function_max_F_M_vol(Optimize_theta, Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

%% Initial test to verify the consistence of the input:
size_beta = size(beta);
size_theta = size(theta);
if max(size_beta) ~= n &&  max(size_theta) ~= n
    fprintf('Arm angles defined not consistent with the number of arms')
    return;
end

%% Compute different metrix for a drone design
step = 0.25;
[~, ~, ~, ~, ~, ~,~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, F_vol, ~, M_vol]  = Mav_compute_metrics(dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%% Objective function:
% fun = -F_vol -M_vol;
fun = -abs(F_vol) -abs(M_vol);
end
%% Maximize Mmin and Fmin and minimize the inertia
function [fun] = objective_function_max_M_F_in_xyz(Optimize_theta, Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end

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
[m, ~, ~, ~, ~, ~] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
D =  [0 0 1; 0 1 0; 1 0 0];
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
            if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Fdes = 0*d;
               break;
            end
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
        else % If alpha0 and w0 does not respect their bounds anymore.
            if i == 1 % if there are no viable solution 
               w0 = 0*w0;
               alpha0 = 0*alpha0;
               Mdes = 0*d;
               break;
            end
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

%% Objecticve function
fun = -sum(Mmin) -sum(Fmin);
end

%% Nonlinear constraint to force the drone to hover in every direction.
function [c,ceq] = nonlinconF(Optimize_theta,Optimize_L, dec, n, x, L, kf, km, wmin, wmax, alphamin, alphamax, g, max_iterations, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
beta = x(1:n);
if Optimize_theta 
    theta = x(n+1:2*n);
    if Optimize_L 
        L = x(2*n+1);
    end
else
    theta = zeros(1,n);
    if Optimize_L
        L = x(n+1);
    end
end
[m, ~] = Mav_inertias(n, L, theta, beta);
%% Initial test to verify the consistence of the input:
size_beta = size(beta);
size_theta = size(theta);
if max(size_beta) ~= n &&  max(size_theta) ~= n
    fprintf('Arm angles defined not consistent with the number of arms')
    return;
end

%% Compute different metrix for a drone design
step = 0.25;
[~, ~, ~, ~, ~, F] = Mav_compute_metrics(dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%% Objective function:
% fun = -F_vol -M_vol;
c = m*g*ones(size(vecnorm(F))) -vecnorm(F);
ceq = [];
end

end
