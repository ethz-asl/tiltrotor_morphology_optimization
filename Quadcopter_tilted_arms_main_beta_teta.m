
%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
clear all;
close all;
%% Design parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
g = 9.81;
m = Mb+4*Mp;% drone mass [kg]
Ndecimals = 5;
dec = 10.^Ndecimals;
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
alphamin = -pi; 
alphamax = pi; 
alphadotmax = pi;
nmin = 0; % minimum rotor speed allowed [round/s]
nmax =150; % maximum rotor speed allowed [round/s]
betamin = -pi/4; 
betamax = pi/4; 
thetamin = -pi/4; 
thetamax = pi/4; 
step = .25; % 0.1, 0.2, 0.25, 0.5, 1
beta0 = [0 0 0 0];
theta0 = [0 0 0 0];

%Parameters for the optimization of alpha and n:
opt_iterations = 150; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true;
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqp-legacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;
beta = beta0;
theta = theta0;
for i = 2:opt_iterations % loop that performs the optimization until the solution is the best possible.
    [beta(i, :), theta(i, :), exitflag] = Quadcopter_tilted_arms_optimize_beta_teta(kf, km, L, g, Mb, Mp, R, nmin, nmax, betamin, betamax, thetamin, thetamax, alphadotmax, alphamin, alphamax, opt_iterations, step, beta(i-1,:), theta(i-1,:), Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);
    if isequal(round(beta(i,:)*10^4)/10^4,round(beta(i-1,:)*10^4)/10^4) && isequal(round(theta(i,:)*10^4)/10^4,round(theta(i-1,:)*10^4)/10^4)
        break;
    end
end
beta = beta(end,:);
theta = theta(end,:);
[m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, eye(3), [0 0 0 0], beta, theta, [0 0 0 0].', L, g, Mb, Mp, R, false);
[D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, opt_iterations, alphadotmax);
Quadcopter_tilted_arms_plot(1, 1, theta, beta,  D, F, Feff, M,Meff, Heff, C, L, R, Op1, Op2, Op3, Op4, step, worthF, worthM, worthH, worthC, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
A1 = [1, rad2deg(beta(1)), rad2deg(beta(2)), rad2deg(beta(3)), rad2deg(beta(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax, Hmin, Hmax];
formatSpec = 'Design %3.0f with optimization on α (β = [%2.0f %2.0f %2.0f %2.0f], θ = [%2.0f %2.0f %2.0f %2.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f, Hmin = %2.2f Hmax = %2.2f\n';
fprintf(formatSpec, A1);