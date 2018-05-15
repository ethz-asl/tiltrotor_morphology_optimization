clear all;
close all;
%% Design parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
g = 9.81; 
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
wmin = 0; % minimum rotor speed allowed [round/s]
wmax =150; % maximum rotor speed allowed [round/s]
alphamin = -pi; % minimum  tilting angle allowed [rad]
alphamax = pi; % maximum tilting angle allowed [rad]
alphadotmax = pi; % max speed of the rotor tilting [rad/s]

step = .1; % 0.1, 0.2, 0.25, 0.5, 1

dec = 10^5; % decimal to round the values returned by the simulation

%Parameters for the optimization of alpha and n:
max_iterations = 200; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true;
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqplegacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;

% [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(4, kf, km, eye(3), [0 0 0 0], [pi/6, pi/6, pi/6, pi/6], [0 0 0 0], [0 0 0 0].', L, g, Mb, Mp, R, dec, false);
% [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, 4, [pi/6, pi/6, pi/6, pi/6] ,[0 0 0 0], L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
% Mav_plot(4, wRb, 25, 25, [0 0 0 0], [pi/6, pi/6, pi/6, pi/6],  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
% 
% [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(4, kf, km, eye(3), [0 0 0 0], [-pi/6, -pi/6, -pi/6, -pi/6], [0 0 0 0], [0 0 0 0].', L, g, Mb, Mp, R, dec, false);
% [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, 4, [-pi/6, -pi/6, -pi/6, -pi/6] ,[0 0 0 0], L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
% Mav_plot(4, wRb, 26, 26, [0 0 0 0], [-pi/6, -pi/6, -pi/6, -pi/6],  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)


%% Test to evaluate simulation with a variable number of propeller n:
beta = [0 0];
A = [];
Formatspecs = [];

parfor i = 3:7
    tStart = tic;
    
    A1 = i;
    formatSpec = 'Beginning optimizatin for design %d\nComputing...\n';
    fprintf(formatSpec, A1);
    
    beta = zeros(1,n);
    theta = beta;
    [rows, column] = size(beta); 
    n = max([rows, column]); % Number of propellers
    
    Mb = 0.25*n; % mass body [kg]
    R = 0.01*n; % Radius of the body (Body assumed to be a sphere)
         
    [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, Mb, Mp, R, dec, false);
    [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, 2*i-1, i, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, 2*i, i, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    
    tEnd = toc(tStart);
    A1= [i, floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Finished optimizing for design %d in %d minutes and %2.2f seconds\n';
    fprintf(formatSpec, A1);
end

% %% Test to evaluate the need of optimization for α and n:
% beta = pi/6*[0 0 0 0; 1 -1 -1 1; -0.5 -0.5 0.5 0.5; 1 1 1 1; -0.5 -0.5 -0.5 -0.5; 0 0 0 1; 1 0 1 0; 0 0.5 0.5 0; 1 0 -1 1].';
% A = [];
% Formatspecs = [];
% for i = beta
%     theta = [0 0 0 0]; %[pi/18 pi/18 pi/18 pi/18];
%     [rows, column] = size(i); % Number of propellers
%     n = max([rows, column]);
%     A1 = ii;
%     formatSpec = 'Beginning optimizatin for design %3.0f\nComputing...\n';
%     fprintf(formatSpec, A1);
%     [m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), [0 0 0 0], i, theta, [0 0 0 0].', L, g, Mb, Mp, R, dec, false);
%     [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, i ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%     Mav_plot(n, wRb, 2*ii-1, ii, theta, i,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
%     A1 = [ii, rad2deg(i(1)), rad2deg(i(2)), rad2deg(i(3)), rad2deg(i(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax, Hmin, Hmax];
%     A = [A; A1];
%     formatSpec = 'Design %3.0f without optimization on α (β = [%2.0f %2.0f %2.0f %2.0f], θ = [%2.0f %2.0f %2.0f %2.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f, Hmin = %2.2f Hmax = %2.2f\n';
%     Formatspecs = [Formatspecs; formatSpec];
%     [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, i ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%     Mav_plot(n, wRb, 2*ii, ii, theta, i,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
%     A1 = [ii, rad2deg(i(1)), rad2deg(i(2)), rad2deg(i(3)), rad2deg(i(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax, Hmin, Hmax];
%     A = [A; A1];
%     formatSpec = 'Design %3.0f with optimization on α    (β = [%2.0f %2.0f %2.0f %2.0f], θ = [%2.0f %2.0f %2.0f %2.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f, Hmin = %2.2f Hmax = %2.2f\n';
%     Formatspecs = [Formatspecs; formatSpec];
%     A1 = ii;
%     formatSpec = 'Finished optimizing for design %3.0f\n';
%     fprintf(formatSpec, A1);
%     ii = ii + 1;
% end
% for i = 1:ii-1
%     fprintf(Formatspecs(i, :), A(i,:));
% end