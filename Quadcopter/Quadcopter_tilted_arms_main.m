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
nmin = 0; % minimum rotor speed allowed [round/s]
nmax =150; % maximum rotor speed allowed [round/s]
alphamin = -pi; % minimum  tilting angle allowed [rad]
alphamax = pi; % maximum tilting angle allowed [rad]
alphadotmax = pi; % max speed of the rotor tilting [rad/s]
ii = 1;
step = .5; % 0.1, 0.2, 0.25, 0.5, 1
theta = [0 0 0 0];

%Parameters for the optimization of alpha and n:
opt_iterations = 200; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true;
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqp-legacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;
%% Find design properties:
% beta_step = 1;
% Data = [];
% for i = -1:beta_step:1
%     for j = -1:beta_step:1
%         for k = -1:beta_step:1
%             for l = -1:beta_step:1
%                 if ii ==41 || ii== 25 || ii == 57 || ii ==75 || ii== 73 || ii == 52 || ii == 30 || ii ==9 || ii== 7
%                     plot= true;
%                 else 
%                     plot = false;
%                 end
%                 beta = pi/6*[i j k l];
%                 theta = Quadcopter_tilted_arms_find_theta(beta, optimize_theta);
%                 [Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot, optimize_alpha);
%                 A1 = [ii, rad2deg(beta(1)), rad2deg(beta(2)), rad2deg(beta(3)), rad2deg(beta(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax];
%                 formatSpec = 'Design %3.0f (β = [%3.0f %3.0f %3.0f %3.0f], θ = [%3.0f %3.0f %3.0f %3.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f\n';
%                 fprintf(formatSpec,A1)
%                 ii = ii + 1;
%                 Data = [Data, [Fmin;Fmax;Mmin;Mmax;beta(1);beta(2);beta(3);beta(4)]];
%             end
%         end
%     end
% end

%% Test to evaluate the need of optimization for α and n:
beta = pi/6*[0 0 0 0; 1 -1 -1 1; -0.5 -0.5 0.5 0.5; 1 1 1 1; -0.5 -0.5 -0.5 -0.5; 0 0 0 1; 1 0 1 0; 0 0.5 0.5 0; 1 0 -1 1].';
A = [];
Formatspecs = [];
for i = beta
    theta = [0 0 0 0];
    A1 = ii;
    formatSpec = 'Beginning optimizatin for design %3.0f\nComputing...\n';
    fprintf(formatSpec, A1);
    [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, eye(3), [0 0 0 0], i, theta, [0 0 0 0].', L, g, Mb, Mp, R, false);
%     [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(i ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance,  opt_iterations, alphadotmax);
%     Quadcopter_tilted_arms_plot(wRb, 2*ii-1, ii, theta, i,  D, F, Feff, M,Meff, Heff, C, L, R, Op1, Op2, Op3, Op4, step, worthF, worthM, worthH, worthC, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
%     A1 = [ii, rad2deg(i(1)), rad2deg(i(2)), rad2deg(i(3)), rad2deg(i(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax, Hmin, Hmax];
%     A = [A; A1];
%     formatSpec = 'Design %3.0f without optimization on α (β = [%2.0f %2.0f %2.0f %2.0f], θ = [%2.0f %2.0f %2.0f %2.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f, Hmin = %2.2f Hmax = %2.2f\n';
%     Formatspecs = [Formatspecs; formatSpec];
    [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(i ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance,  opt_iterations, alphadotmax);
    Quadcopter_tilted_arms_plot(wRb, 2*ii, ii, theta, i,  D, F, Feff, M,Meff, Heff, C, L, R, Op1, Op2, Op3, Op4, step, worthF, worthM, worthH, worthC, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    A1 = [ii, rad2deg(i(1)), rad2deg(i(2)), rad2deg(i(3)), rad2deg(i(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax, Hmin, Hmax];
    A = [A; A1];
    formatSpec = 'Design %3.0f with optimization on α    (β = [%2.0f %2.0f %2.0f %2.0f], θ = [%2.0f %2.0f %2.0f %2.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f, Hmin = %2.2f Hmax = %2.2f\n';
    Formatspecs = [Formatspecs; formatSpec];
    A1 = ii;
    formatSpec = 'Finished optimizing for design %3.0f\n';
    fprintf(formatSpec, A1);
    ii = ii + 1;
end
for i = 1:ii-1
    fprintf(Formatspecs(i, :), A(i,:));
end
