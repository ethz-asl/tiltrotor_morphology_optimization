%%%%%%%%%%%% n-copter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
clear all;
close all;
[file_path] = fileparts(mfilename('fullpath'));
addpath(file_path);
file_path = erase(file_path, 'Functions_for_Euler_supercomputer');
addpath([file_path '/Mav_optimization_tool_functions/']);
%% Design parameters
g = 9.81;
Ndecimals = 5;
dec = 10.^Ndecimals;
kf = 3.86e-4; % Propeller thrust coefficient % false[kg.m]
km = 1.5e-5;% Propeller drag coefficient
Lmin = 0.1;
Lmax = 0.5; % Arm length [m]
alphamin = -pi; 
alphamax = pi;
alphadotmax = pi;
wmin = 0; % minimum rotor speed allowed [round/s]
wmax =150; % maximum rotor speed allowed [round/s]
betamin = -4*pi/20;
betamax = 4*pi/10;
thetamin = -pi;
thetamax = pi;
nmin = 3;
nmax = 8;
%% Parameters for the optimization of tilting angles (alpha) and rotor speeds (w):
step = .1; % step to define the number of directions in which to compute forcetorque/hover eff
           % (0.5 -> 98 directions, 0.25 -> 578 directions, 0.1 -> 7490 directions)
max_iterations = 150; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true; % If true performs an optimization on the tilting angles and rotor speeds to max the force/torque/hover eff in every direction
                       % if false uses the angles returned by the static matrix solution
                       
%% Parameters for fmincom fct
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqp-legacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;
beta1 = acos(sqrt(2/3));
%% Flat design for comparison
for n= 5:2:7 % for a n-rotor MAV
    tStart = tic; % start timer
    A1 = n;
    formatSpec = 'Beginning design optimization for a %d-MAV \nComputing...\n';
    fprintf(formatSpec, A1);
    thetamin = -2*pi/n + pi/(4*n);
    thetamax = 2*pi/n - pi/(4*n);
    % Initial solution
    beta = zeros(1,n); % vertical angle of the arms of the n-rotor MAV
    theta = zeros(1,n);
    if mod(n,2) == 0
        for ll = 1:n
            if mod(ll,2) == 0
                beta(ll) = -beta1;
            else
                beta(ll) = beta1;
            end
        end
    else
        beta = ones(1,n)*beta1;
    end
    % L = Lmin + (Lmax-Lmin)/2; % Arm length
    L = Lmax;
    R = n*0.1/4; % Radius of the body (Body assumed to be a sphere)
    
    % Compute metrix for the solution and plot the result
    [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
    [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, optimize_alpha, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, n, n, theta, beta,  D, F, Feff, M,Meff, Heff, L, Op, bRp, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    filename = ['n=' num2str(n) '_opt.fig'];
    saveas(figure(n),filename);
    close(figure(n));
    tEnd = toc(tStart); % log exec time
    
    A1 = [n, floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Design optimization for a %d-MAV finished in %d minutes and %2.2f seconds \n';
    fprintf(formatSpec, A1);
    fprintf(['β = ' mat2str(round(rad2deg(beta)*10^2)/10^2) ', θ = ' mat2str(round(rad2deg(theta)*10^2)/10^2) ', L = ' num2str(L) '\n']);
end