%%%%%%%%%%%% n-copter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
clear all;
close all;
%% Design parameters
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
betamin = -pi/2;
betamax = pi/2;
thetamin = -2*pi/9;
thetamax = 2*pi/9;
nmin = 3;
nmax = 8;
%Parameters for the optimization of alpha and n:
step = .5; % step between points of the force, torque and hover space
max_iterations = 150; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = true;
Algorithm = 'sqp'; %,'sqp' (best tested), 'sqp-legacy' 'interior-point' (way too long), 'active-set'
Display = 'off'; % 'off', 'notify'
maxIter = 10000;
StepTolerance = 1.0000e-6;
ConstraintTolerance = 1.0000e-6;

%% optimize beta, theta and L using the static matrix and fmincom
for n= 3:7
    tStart = tic; % start timer
    A1 = n;
    formatSpec = 'Beginning design optimization for a %d-MAV \nComputing...\n';
    fprintf(formatSpec, A1);
    beta = zeros(1,n);
    for ll = 3:n
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
    R = n*0.1/4; % Radius of the body (Body assumed to be a sphere)
    
%     [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), Lmax, g, dec, false);
%     [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
%     Mav_plot(n, wRb, 2*n-1, 2*n-1, theta, beta,  D, F, Feff, M,Meff, Heff, Lmax, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    
    L = Lmax;
    exitflag = [];
    obj_fun = [];
    out = false;
    for i = 2:max_iterations % loop that performs the optimization until the solution is the best possible.
        [beta(i, :), theta(i, :), L(i), obj_fun(i), exitflag(i)] = Mav_optimize_beta_theta_L(dec, n, kf, km, Lmin, Lmax, L(i-1), g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, beta(i-1, :), theta(i-1, :), Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);                                                                    
        for ii = 1:i
            if isequal(round(beta(i,:)*10^3)/10^3,round(beta(ii,:)*10^3)/10^3) && isequal(round(theta(i,:)*10^3)/10^3,round(theta(ii,:)*10^3)/10^3) && isequal(round(L(i)*10^3)/10^3,round(L(ii)*10^3)/10^3) && i ~= ii
                out = true;
                if obj_fun(i) > obj_fun(ii) && i ~= ii
                    beta(i, :) = beta(ii, :);
                    theta(i, :)= theta(ii, :);
                end
            end
        end
        if out == true
            break;
        end
    end
    beta = round(beta(end,:)*dec)/dec;
    theta = round(theta(end,:)*dec)/dec;
    L = round(L(end)*dec)/dec;
    [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
    [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, n, n, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    tEnd = toc(tStart);
    A1 = [n, floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Design optimization for a %d-MAV finished in %d minutes and %2.2f seconds \n';
    fprintf(formatSpec, A1);
    fprintf(['β = ' mat2str(round(rad2deg(beta)*10^2)/10^2) ', θ = ' mat2str(round(rad2deg(theta)*10^2)/10^2) ', L = ' num2str(L) '\n']);
end

% %% optimize n, beta, theta and L using ga
% tStart = tic; % start timer
% fprintf('Beginning design optimization for a n-MAV \nComputing...\n');
% [n, beta, theta, L, obj_fun, exitflag] = Mav_optimize_n_beta_theta_L(dec, kf, km, nmin, nmax, Lmin, Lmax, g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, Display, ConstraintTolerance);
% beta = round(beta(end,:)*dec)/dec;
% theta = round(theta(end,:)*dec)/dec;
% L = round(L(end)*dec)/dec;
% R = n*0.1/4; % Radius of the body (Body assumed to be a sphere)
% [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
% [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, true, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
% Mav_plot(n, wRb, n, n, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
% tEnd = toc(tStart);
% A1 = [floor(tEnd/60), rem(tEnd,60)];
% formatSpec = 'Design optimization for a n-MAV finished in %d minutes and %2.2f seconds \n';
% fprintf(formatSpec, A1);