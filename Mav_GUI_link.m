function [] = Mav_GUI_link(cost_fct_case, Optimize_theta, Optimize_L, Optimize_n, direction, design_number, nmin, nmax, Lmin, Lmax, betamin, betamax, thetamin, thetamax, n, L, theta, beta, max_iterations, step, optimize_alpha, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)

%% Initialize some parameters
[g, dec, kf, km, alphamin, alphamax, wmin, wmax] = Mav_parameters();
closereq
if Optimize_n
    %% optimize number of propellers (n), arms vertical angles (beta), horizontal angles (theta) and length (L)
    tStart = tic; % start timer
    fprintf('Beginning design optimization for a n-MAV \nComputing...\n');
    [n, beta, theta, L, obj_fun, exitflag] = Mav_optimize_n_beta_theta_L(cost_fct_case, Optimize_theta, Optimize_L, dec, L, direction, kf, km, nmin, nmax, Lmin, Lmax, g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, Display, ConstraintTolerance)
    beta = round(beta(end,:)*dec)/dec;
    theta = round(theta(end,:)*dec)/dec;
    L = round(L(end)*dec)/dec;
    R = n*0.1/4; % Radius of the body (Body assumed to be a sphere)
    [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
    [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, optimize_alpha, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, design_number, design_number, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    tEnd = toc(tStart);
    A1 = [floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Design optimization for a n-MAV finished in %d minutes and %2.2f seconds \n';
    fprintf(formatSpec, A1);
else
    tStart = tic; % start timer
    A1 = n;
    formatSpec = 'Beginning design optimization for a %d-MAV \nComputing...\n';
    fprintf(formatSpec, A1);
    R = n*0.1/4; % Radius of the body (Body assumed to be a sphere)
    % Perform the optimization on the n-rotor MAV design
    exitflag = [];
    obj_fun = [];
    out = false;
    for i = 2:max_iterations % loop that performs the optimization until the solution is the best possible.
        
        % As an initial solution feed fmincom with the solution of the last iteration
        [beta(i, :), theta(i, :), L(i), obj_fun(i), exitflag(i)] = Mav_optimize_beta_theta_L(cost_fct_case, Optimize_theta, Optimize_L, direction, dec, n, kf, km, Lmin, Lmax, L(i-1), g, wmin, wmax, betamin, betamax, thetamin, thetamax, alphamin, alphamax, max_iterations, beta(i-1,:), theta(i-1,:), Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);

        % Loop to test if the found solution in this iteration is the same
        % as the one in previous iteration. If so -> converged -> break initial loop
        for ii = 1:i
            if isequal(round(beta(i,:)*10^3)/10^3,round(beta(ii,:)*10^3)/10^3) && i ~= ii % && isequal(round(theta(i,:)*10^3)/10^3,round(theta(ii,:)*10^3)/10^3) && isequal(round(L(i)*10^3)/10^3,round(L(ii)*10^3)/10^3)
                out = true;
                if obj_fun(i) > obj_fun(ii) && i ~= ii
                    beta(i, :) = beta(ii, :);
                    theta(i, :)= theta(ii, :);
                    L(i) = L(i-1);
                end
            end
        end
        if out == true
            break;
        end
    end
    
    % Solution when algorithm has converged:
    beta = round(beta(end,:)*dec)/dec;
    theta = round(theta(end,:)*dec)/dec;
    L = round(L(end)*dec)/dec;
    
    % Compute metrix for the solution and plot the result
    [~, ~, ~, ~, Op, bRp] = Mav_dynamic(n, kf, km, eye(3), zeros(n,1), beta, theta, zeros(n,1), L, g, dec, false);
    [wRb, D, Heff, ~, ~, F,~, ~, Feff, M, ~, ~, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, Lmax, kf, km, wmin, wmax, alphamin, alphamax, g, step, optimize_alpha, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Mav_plot(n, wRb, design_number, design_number, theta, beta,  D, F, Feff, M,Meff, Heff, L, R, Op, bRp, step, worthF, worthM, worthH, number_of_directions, true, TRI, F_surf, F_vol, M_surf, M_vol)
    
    tEnd = toc(tStart); % log exec time
    
    A1 = [n, floor(tEnd/60), rem(tEnd,60)];
    formatSpec = 'Design optimization for a %d-MAV finished in %d minutes and %2.2f seconds \n';
    fprintf(formatSpec, A1);
    fprintf(['β = ' mat2str(round(rad2deg(beta)*10^2)/10^2) ', θ = ' mat2str(round(rad2deg(theta)*10^2)/10^2) ', L = ' num2str(L) '\n']);

end
end