function [betastar, tetastar, Lstar, exitflag] = Mav_optimize_beta_theta(dec, n, kf, km, Lmin, Lmax, L0, g, Mb, Mp, m, R, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)                                                                        
% [betastar, tetastar, exitflag] = Mav_optimize_beta_theta(n, kf, km, L, g, Mb, Mp, R, nmin, nmax, betamin, betamax, thetamin, thetamax, alphadotmax, alphamin, alphamax, opt_iterations, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%MAV_OPTIMIZE_BETA_THETA find optimal tilting angles and Rotor speed  
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
lb_beta = betamin*ones(n,1);
ub_beta = betamax*ones(n,1);
lb_theta = thetamin*ones(n,1);
ub_theta = thetamax*ones(n,1);
lb_L = Lmin;
ub_L = Lmax;
lb = [lb_beta; lb_theta; lb_L];% lower bound
ub = [ub_beta; ub_theta; ub_L];% upper bound 

% initial guess:
x0 = [beta0 theta0 L0];
         
% optimization options       
options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

% actual optimization
[xstar,~,exitflag,~] = fmincon(@ (x) objective_function(dec, n, x, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations), x0, A, b, Aeq, beq, lb, ub,[],  options);

% Solution of the optimization
betastar = xstar(1:n);
tetastar = xstar(n+1:2*n);
Lstar = xstar(2*n+1);

%% Objective function function  
function [fun] = objective_function(dec, n, x, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations)
    beta = x(1:n);
    theta = x(n+1:2*n);
    L = x(2*n+1);
    [~, ~, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, ~, ~, ~, ~, ~, ~, F_vol, ~, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    fun = -Hmin - Fmin - Mmin;
end
end
