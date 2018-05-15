%%%%%%%%%%%% n-copter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
clear all;
close all;
%% Design parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
g = 9.81;
Ndecimals = 5;
dec = 10.^Ndecimals;
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
Lmin = 0.1;
Lmax = 0.5; % Arm length [m]
alphamin = -pi; 
alphamax = pi; 
alphadotmax = pi;
wmin = 0; % minimum rotor speed allowed [round/s]
wmax =150; % maximum rotor speed allowed [round/s]
betamin = -pi/4; 
betamax = pi/4; 
thetamin = -2*pi/9; 
thetamax = 2*pi/9;
nmin = 3;
nmax = 12;
%Parameters for the optimization of alpha and n:
step = .25;
max_iterations = 150; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true;
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqp-legacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;

% [m, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta0, theta0, zeros(n,1), Lmax, g, Mb, Mp, R, dec, false);
% [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta0 ,theta0, Lmax, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
% Mav_plot(n, wRb, 1, 1, theta0, beta0,  D, F, Feff, M,Meff, Heff, Lmax, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)

%% optimize beta and theta using the static matrix and fmincom
for n= 7:-1:3
    step = .25;
    tStart = tic;
    A1 = n;
    formatSpec = 'Beginning design optimization for a %d-MAV \nComputing...\n';
    fprintf(formatSpec, A1);
    beta = zeros(1,n);
    for ll = 1:n
        if mod(ll,2) == 0
            beta(ll) = -pi/5;
        else
            beta(ll) = pi/5;
        end
    end
    if n==3
        beta = ones(1,n)*pi/5;
    end
    theta = zeros(1,n);

%     beta = zeros(1,n);
%     theta = zeros(1,n);
    L = 0.15;
    exitflag = [];
    m = n*Mp + Mb;
    for i = 2:max_iterations % loop that performs the optimization until the solution is the best possible.
        [beta(i, :), theta(i, :), L(i), exitflag(i)] = Mav_optimize_beta_theta(dec, n, kf, km, Lmin, Lmax, L(i-1), g, Mb, Mp, m, R, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, step, beta(i-1, :), theta(i-1, :), Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);                                                                    
        if isequal(round(beta(i,:)*10^4)/10^4,round(beta(i-1,:)*10^4)/10^4) && isequal(round(theta(i,:)*10^4)/10^4,round(theta(i-1,:)*10^4)/10^4) && isequal(round(L(i)*10^4)/10^4,round(L(i-1)*10^4)/10^4)
            break;
        end
    end
    beta = round(beta(end,:)*dec)/dec;
    theta = round(theta(end,:)*dec)/dec;
    L = round(L(end)*dec)/dec;
    step = .1;
    [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, Mb, Mp, R, dec, false);
    [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, 2*n, 2*n, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    tEnd = toc(tStart);
    A1 = [n, floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Design optimization for a %d-MAV finished in %d minutes and %2.2f seconds \n';
    fprintf(formatSpec, A1);
    fprintf(['β = ' mat2str(beta) ', θ = ' mat2str(theta) ', L = ' num2str(L) '\n']);
end