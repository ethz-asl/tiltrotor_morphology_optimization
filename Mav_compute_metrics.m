function [wRb, D_unit2, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, length_D, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, step, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations)
%MAV_COMPUTE_METRICS computes a lot of metrics for a given design of Mav
%   Design defined by the number of arms and their angles (beta & theta) and other parameters

%%%%%%%%%%%% MAV with tilting rotor and tilted arms design optimization%%%%%%%%%%%%

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
%
%         ...
%
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
%% Loop to create the direction matrix D with the desired number of direction:
D = zeros(3,(2/step+1)^3); % Matrix containing every direction for which we want to compute the metrics
length_D = 0;
for i = 1:-step:-1
    for j = 1:-step:-1
        for k = 1:-step:-1% number of directions depends only on step size
            length_D = length_D+1;
            D(:,length_D) = [i j k].';
        end
    end
end
i0 = ~vecnorm(D); 
D(:,i0) = [];% Eliminate the direction [0; 0; 0] 
D_unit = D./vecnorm(D); % Create a normalized matrix of direction.
D_unit = round(D_unit*dec)/dec; % Round D_unit
[D_unit,ia,~] = unique(D_unit.', 'stable', 'rows'); % Eliminate redundant directions in normalized D
D_unit = D_unit.';
D = D(:,ia.'); % Eliminate redundant directions in D
%% initialization of the parameters:
[~, length_D] = size(D);

F = []; % Matrix containing the maximum force appliable by the design in every direction of D
Feff = zeros(1, length_D); % Vector containing the efficiency of the drone when applying the maximum force in every direction of D
alphastarF = [];% Vector containing the optimal tilting angles to obtain the max force in every direction of D
wstarF = []; % Vector containing the optimal propeller speed to obtain the max force in every direction of D

M = []; % Matrix containing the maximum torques appliable by the design in every direction of D
Meff = zeros(1, length_D);% Vector containing the efficiency of the drone when applying the maximum torque in every direction of D
alphastarM = [];% Vector containing the optimal tilting angles to obtain the max torque in every direction of D
wstarM = []; % Vector containing the optimal propeller speed to obtain the max torque in every direction of D

Heff = []; % Vector containing the efficiency of the drone when hovering with the weight oriented in direction -D 
alphastarH = [];% Vector containing the optimal tilting angles to hover in ev ery direction of D
wstarH = []; % Vector containing the optimal propeller speed to obtain the max force in every direction of D

D_unit2 = [];

worthF = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal force
worthM = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal torque
worthH = 0; % Counter to quantify the efficiency of the fmincom optimization on hover efficiency

%% Test parameters:
iiif = 0; % incremented if the optimization returns a worse solution than the initial to the maximum force problem in one direction
iiim = 0; % incremented if the optimization returns a worse solution than the initial to the maximum torque problem in one direction
iiih = 0; % incremented if the optimization returns a worse solution than the initial to the optimal hover problem in one direction

% exitflagf = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the maximal force in direction d
% exitflagm = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the maximal torque in direction d
% exitflagh = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the optimal hover in direction d

% successh = 0;
% successf = 0;
% successm = 0;
% failh = 0;
% failf = 0;
% failm = 0;

[m, ~] = Mav_inertias(n, L, theta, beta);

%% Loop to perform the optimizations of the max force, torque, hover efficiency in every direction
for ii = 1:1:length_D
    d = D_unit(:,ii);
    D_unit2(:,ii) = d;
    %% First, find the max force in direction d using static matrix
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
    
    % calculate linear acceleration with this alphastar and nstar
    [m, ~, pdotdot] = Mav_dynamic(n, kf, km, wRb, alpha0, beta, theta,w0, L, g, dec, false);
    F0 = m*pdotdot;
    FN0 = norm(F0);
%     exitflag1 = 0;
    %% find the max force in direction d using fmincom and static matrix solution as initial solution
    if optim % performs the optimisation only if optim is true                
        alphastar = [alpha0 alpha0];
        wstar = [w0 w0];
        Fstar = [F0 F0];
        FNstar = [FN0 FN0];
        i = 3;
        % loop that performs the optimization until the solution is the
        % best possible feeding fmincom with the previous solution as starting point
        while( i < max_iterations) 
            if FNstar(i-1) ~= 0
                % Perform the optimization and find the max thrust in direction d
                [alphastarloop, nstarloop, ~] = Mav_maximize_force(n, kf, wmin, wmax, alphamin, alphamax, alphastar(:,i-1), wstar(:,i-1), d, beta, theta, wRb, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            else
                % Perform the optimization and find the max thrust in direction d
                [alphastarloop, nstarloop, ~] = Mav_maximize_force(n, kf, wmin, wmax, alphamin, alphamax, alphastar(:,i-2), wstar(:,i-2), d, beta, theta, wRb, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            end
            wstar(:,i) = round(dec*nstarloop)/dec; % rotors speeds after optimization
            alphastar(:, i) = round(dec*alphastarloop)/dec; % rotors orientations after optimization

            % Calculate angular and linear acceleration with this alphastar and nstar
            [m, ~, pdotdot]  = Mav_dynamic(n, kf, km, wRb, alphastar(:, i), beta, theta, wstar(:,i), L, g, dec, false);
            pdotdot = round(dec*pdotdot)/dec;
            Fstar(:,i) = m*pdotdot; % Force applied to the body with the propellers in this 
            FNstar(i) = norm(Fstar(:,i));
            
            % If this solution breaks the constraint: Fstar parallel to d
            % Return to the previous solution and quit the loop
            if ~isequal(round((Fstar(:,i)/FNstar(i))*(dec/1000))/(dec/1000),round(d*(dec/1000))/(dec/1000))
                alphastar(:,i) = alphastar(:,i-2);
                wstar(:,i) = wstar(:,i-2);
                Fstar(:,i) = Fstar(:,i-2);
                FNstar(i) = FNstar(i-2);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j,i)<alphamin || alphastar(j,i)>alphamax
                    alphastar(:,i) = alphastar(:, i-2);
                    wstar(:,i) = wstar(:,i-2);
                    Fstar(:,i) = Fstar(:,i-2);
                    FNstar(i) = FNstar(i-2);
                    break;
                end
            end
            
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(FNstar(i)*dec)/dec < round(FNstar(i-2)*dec)/dec
                alphastar(:, i) = alphastar(:, i-2);
                wstar(:,i) = wstar(:,i-2);
                Fstar(:,i) = Fstar(:,i-2);
                FNstar(i) = FNstar(i-2);
                Fdes = d*n*wmax^2*kf;
                Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
                % Retrieve rotors speeds and orientations from Fdec
                [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec);
                % Test to see if this wi, alph have already
                test = [];
                [~, columns] = size(alphastar);
                for j = 1:columns
                    if isequal(round(alphai*dec)/dec,round(alphastar(:,j)*dec)/dec) && isequal(round(wi*dec)/dec, round(wstar(:, j)*dec)/dec)
                        test = [test true];
                    end
                end
                if isempty(test) && i < max_iterations-2
                    i = i+1;
                    alphastar(:,i) = alphai;
                    wstar(:,i) = wi;
                    Fstar(:,i) = Fstar(:,i-1);
                    FNstar(i) = FNstar(i-1);
                    i = i+1;
                    continue;
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarF, wstarF, alphastar, wstar, FNstar(i), vecnorm(F), dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Fstar(:,i) = Fstar(:,i-1);
                        FNstar(i) = FNstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
%             exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(FNstar(i)*(dec/100))/(dec/100) == round(FNstar(i-2)*(dec/100))/(dec/100)
                if round(FNstar(i)*dec)/dec <= round(FN0*dec)/dec
                    Fdes = d*n*wmax^2*kf;
                    Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
                    % Retrieve rotors speeds and orientations from Fdec
                    [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec);
                    % Test to see if this wi, alph have already
                    test = [];
                    [~, columns] = size(alphastar);
                    for j = 1:columns
                        if isequal(round(alphai*dec)/dec,round(alphastar(:,j)*dec)/dec) && isequal(round(wi*dec)/dec, round(wstar(:, j)*dec)/dec)
                            test = [test true];
                        end
                    end
                    if isempty(test) && i < max_iterations-2
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Fstar(:,i) = Fstar(:,i-1);
                        FNstar(i) = FNstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarF, wstarF, alphastar, wstar, FNstar(i), vecnorm(F), dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Fstar(:,i) = Fstar(:,i-1);
                        FNstar(i) = FNstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
            if i < max_iterations-2
                i = i+1;
                alphastar(:,i) = zeros(size(alpha0));
                wstar(:,i) = zeros(size(w0));
                Fstar(:,i) = [0; 0; 0];
                FNstar(i) = 0;
                i = i+1;
            else
                break;
            end
        end
%         exitflagf(ii) = exitflag1;
        if round(FNstar(end)*(dec/100))/(dec/100) < round(FN0*(dec/100))/(dec/100)
            iiif = iiif +1;
        end
        F(:,ii) = Fstar(:,end);% Force produced by the MAV
        Feff(ii) = kf*(norm(wstar(:,end))^2); % Efficiency of this force
        alphastarF(:,ii) = alphastar(:,end);% rotors orientations needed to obtain that force
        wstarF(:,ii) = wstar(:,end);% rotors speed needed to obtain that force
    else
        % If optim = false -> static matrix sol:
        F(:,ii) = F0;% Force produced by the MAV
        Feff(ii) = kf*(norm(w0)^2); % Efficiency of this force
        alphastarF(:,ii) = alpha0;% rotors orientations needed to obtain that force
        wstarF(:,ii) = w0;% rotors speed needed to obtain that force
    end
    
    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(norm(F(:,ii))*(dec/100))/(dec/100) <= round(FN0*(dec/100))/(dec/100)
        worthF = worthF+1;
    end
    
    %% find max torque in direction d using static matrix
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
    MN0 = norm(M0);
    alphastar = [alpha0 alpha0];
    wstar = [w0 w0];
    Mstar = [M0 M0];
    MNstar = [MN0 MN0];
%     exitflag1 = 0;
    %% find the max torque in direction d using fmincom and static matrix solution as initial solution
    if optim
        i = 3;
        % loop that performs the optimization until the solution is the
        % best possible feeding fmincom with the previous solution as starting point
        while( i < max_iterations)
            if MNstar(i-1) ~= 0
                % Perform the optimization and find the max torque in direction d
                [alphastarloop, nstarloop, ~] = Mav_maximize_torque(n, kf, km, L, wmin, wmax, alphamin, alphamax, alphastar(:,i-1), wstar(:,i-1), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            else
                % Perform the optimization and find the max torque in direction d
                [alphastarloop, nstarloop, ~] = Mav_maximize_torque(n, kf, km, L, wmin, wmax, alphamin, alphamax, alphastar(:,i-2), wstar(:,i-2), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            end
            wstar(:,i) = round(dec*nstarloop)/dec; % optimal rotors speeds
            alphastar(:, i) = round(dec*alphastarloop)/dec; % optimal rotors orientations
            
            % Calculate angular acceleration with this alphastar and nstar
            [~, Ib, ~, wbdot] = Mav_dynamic(n, kf, km, wRb, alphastar(:,i), beta, theta, wstar(:,i), L, g, dec, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar(:,i) = Ib*wbdot;
            MNstar(i) = norm(Mstar(:,i));
            
            % If this solution breaks the constraint: Mstar parallel to d
            % Return to the previous solution and quit the loop
            if ~isequal(round((Mstar(:,i)/MNstar(i))*(dec/1000))/(dec/1000),round(d*(dec/1000))/(dec/1000))
                wstar(:,i) = wstar(:,i-2);
                alphastar(:, i) = alphastar(:, i-2);
                Mstar(:,i) = Mstar(:,i-2);
                MNstar(i) = MNstar(i-2);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j, i)<alphamin || alphastar(j, i)>alphamax
                    wstar(:,i) = wstar(:,i-2);
                    alphastar(:, i) = alphastar(:, i-2);
                    Mstar(:,i) = Mstar(:,i-2);
                    MNstar(i) = MNstar(i-2);
                    break;
                end
            end
            
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(MNstar(i)*dec)/dec < round(MNstar(i-2)*dec)/dec
                alphastar(:, i) = alphastar(:, i-2);
                wstar(:,i) = wstar(:,i-2);
                Mstar(:,i) = Mstar(:,i-2);
                MNstar(i) = MNstar(i-2);
                Mdes = d*L*n*wmax^2*kf;
                Fdec = A_M_staticinv*Mdes; % Fdec = inv(Astatic)*Fdes
                % Retrieve rotors speeds and orientations from Fdec
                [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec);
                % Test to see if this wi, alph have already
                test = [];
                [~, columns] = size(alphastar);
                for j = 1:columns
                    if isequal(round(alphai*dec)/dec,round(alphastar(:,j)*dec)/dec) && isequal(round(wi*dec)/dec, round(wstar(:, j)*dec)/dec)
                        test = [test true];
                    end
                end
                if isempty(test) && i < max_iterations-2
                    i = i+1;
                    alphastar(:,i) = alphai;
                    wstar(:,i) = wi;
                    Mstar(:,i) = Mstar(:,i-1);
                    MNstar(i) = MNstar(i-1);
                    i = i+1;
                    continue;
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarM, wstarM, alphastar, wstar, MNstar(i), vecnorm(M), dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Mstar(:,i) = Mstar(:,i-1);
                        MNstar(i) = MNstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
%             exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(MNstar(i)*(dec/100))/(dec/100) == round(MNstar(i-2)*(dec/100))/(dec/100)
                % if the optimization did not provide a better solution
                % than the static matrix one
                if round(MNstar(i)*dec)/dec <= round(MN0*dec)/dec
                    %Test an overestimated solution (max theoretical solution)
                    Mdes = d*L*n*wmax^2*kf;
                    Fdec = A_M_staticinv*Mdes; % Fdec = inv(Astatic)*Fdes
                    [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec); % Retrieve rotors speeds and orientations from Fdec
                    % Test to see if this wi, alphi have already beentested
                    % (avoid infinite loop)
                    test = [];
                    [~, columns] = size(alphastar);
                    for j = 1:columns
                        if isequal(round(alphai*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(wi*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                            test = [test true];
                        end
                    end
                    % if alphai and wi have not been tested
                    if isempty(test) && i < max_iterations-2
                        i = i+1;
                        alphastar(:,i) = alphai; % replace the last alphastar by alphai
                        wstar(:,i) = wi; % replace the last wstar by wi
                        Mstar(:,i) = Mstar(:,i-1);
                        MNstar(i) = MNstar(i-1);
                        i = i+1;
                        continue; % and retry the optimizattion
                    end
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarM, wstarM, alphastar, wstar, MNstar(i), vecnorm(M), dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Mstar(:,i) = Mstar(:,i-1);
                        MNstar(i) = MNstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
            if i < max_iterations-2
                i = i+1;
                alphastar(:,i) = zeros(size(alpha0));
                wstar(:,i) = zeros(size(w0));
                Mstar(:,i) = [0; 0; 0];
                MNstar(i) = 0;
                i = i+1;
            else
                break;
            end
        end
        % Check if the optimized solution is not worse than the initial
        if round(MNstar(end)*(dec/100))/(dec/100) < round(MN0*(dec/100))/(dec/100)
            iiim = iiim +1;
        end
%         exitflagm(ii) = exitflag1;
    end
    M(:,ii) = Mstar(:,end);% Torque produced by the MAV
    Meff(ii) = L*kf*(norm(wstar(:,end))^2); % Efficiency of this torque 
    alphastarM(:,ii) = alphastar(:, end); % rotors orientations needed to obtain that torque
    wstarM(:,ii) = wstar(:,end);% rotors speed needed to obtain that torque

    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(MNstar(end)*(dec/100))/(dec/100) <= round(MN0*(dec/100))/(dec/100)
        worthM = worthM+1;
    end
    
    %% find hover efficiency in direction d using static matrix
    % find initial alpha and n for the optimisation to find the best hover in orientation d
    Fdes = m*g*d;
    Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
    
    % Retrieve rotors speeds and orientations from Fdec
    [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
    
    H0 = m*g/(kf*norm(w0)^2);
    H0 = round(dec*H0)/dec;
    wstar = [w0 w0];
    alphastar = [alpha0 alpha0];
    Hstar = [H0 H0];
%     exitflag1 = 0;
    %% find the best hover efficiency in direction d using fmincom and static matrix solution as initial solution
    if optim
        i = 3;        
        % loop that performs the optimization until the solution is the
        % best possible feeding fmincom with the previous solution as starting point
        while( i < max_iterations)
            if Hstar(i-1) ~= 0
                % Perform the optimization and find the best hover in direction d
                [alphastarloop, nstarloop, ~, TH] = Mav_optimize_hover(n, kf, Fdes, wmin, wmax, alphamin, alphamax, alphastar(:, i-1), wstar(:,i-1), beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);
            else
                % Perform the optimization and find the best hover in direction d
                [alphastarloop, nstarloop, ~, TH] = Mav_optimize_hover(n, kf, Fdes, wmin, wmax, alphamin, alphamax, alphastar(:, i-2), wstar(:,i-2), beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);
            end
            wstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(:, i) = round(dec*alphastarloop)/dec;
            Hstar(i) = m*g/(kf*norm(wstar(:,i))^2);
            Hstar(i) = round(dec*Hstar(i))/dec;
            
            % If this solution breaks the constraint: TH = Fdes
            % Return to the previous solution and quit the loop
            if ~isequal(round(TH*(dec/1000))/(dec/1000),round(Fdes*(dec/1000))/(dec/1000))
                alphastar(:, i) = alphastar(:, i-2);
                wstar(:,i) = wstar(:,i-2);
                Hstar(i) = Hstar(i-2);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j, i)<alphamin || alphastar(j, i)>alphamax
                    alphastar(:, i) = alphastar(:, i-2);
                    wstar(:,i) = wstar(:,i-2);
                    Hstar(i) = Hstar(i-2);
                    break;
                end
            end            
           
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(Hstar(i)*dec)/dec < round(Hstar(i-2)*dec)/dec
                alphastar(:, i) = alphastar(:, i-2);
                wstar(:,i) = wstar(:,i-2);
                Hstar(i) = Hstar(i-2);
                %Test an underestimate solution (max theoretical solution)
                Fdes2 = d*n*(wmax-wmin)^2/400*kf;
                Fdec = A_F_staticinv*Fdes2; % Fdec = inv(Astatic)*Fdes
                % Retrieve rotors speeds and orientations from Fdec
                [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec);
                % Test to see if this wi, alph have already
                test = [];
                [~, columns] = size(alphastar);
                for j = 1:columns
                    if isequal(round(alphai*dec)/dec,round(alphastar(:,j)*dec)/dec) && isequal(round(wi*dec)/dec, round(wstar(:, j)*dec)/dec)
                        test = [test true];
                    end
                end
                if isempty(test) && i < max_iterations-2
                    i = i+1;
                    alphastar(:,i) = alphai;
                    wstar(:,i) = wi;
                    Hstar(i) = Hstar(i-1);
                    i = i+1;
                    continue;
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarH, wstarH, alphastar, wstar, Hstar(i), Heff, dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Hstar(i) = Hstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
%             exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(Hstar(i)*(dec/100))/(dec/100) == round(Hstar(i-2)*(dec/100))/(dec/100)
                if round(Hstar(i)*dec)/dec <= round(H0*dec)/dec
                    %Test an underestimate solution (max theoretical solution)
                    Fdes2 = d*n*(wmax-wmin)^2/100*kf;
                    Fdec = A_F_staticinv*Fdes2; % Fdec = inv(Astatic)*Fdes
                    % Retrieve rotors speeds and orientations from Fdec
                    [wi,alphai] = Mav_get_decomposition(n, dec, kf, Fdec);
                    % Test to see if this wi, alph have already
                    test = [];
                    [~, columns] = size(alphastar);
                    for j = 1:columns
                        if isequal(round(alphai*dec)/dec,round(alphastar(:,j)*dec)/dec) && isequal(round(wi*dec)/dec, round(wstar(:, j)*dec)/dec)
                            test = [test true];
                        end
                    end
                    if isempty(test) && i < max_iterations-2
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Hstar(i) = Hstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                if i < max_iterations-2
                    [alphai,wi] = Mav_mirror_points(n, ii, D_unit2, d, alphastarH, wstarH, alphastar, wstar, Hstar(i), Heff, dec);
                    if ~isempty(alphai) && ~isempty(wi)
                        i = i+1;
                        alphastar(:,i) = alphai;
                        wstar(:,i) = wi;
                        Hstar(i) = Hstar(i-1);
                        i = i+1;
                        continue;
                    end
                end
                break;
            end
            if i < max_iterations-2
                i = i + 1;
                alphastar(:, i) = zeros(size(alpha0));
                wstar(:,i) = zeros(size(w0));
                Hstar(i) = 0;
                i = i+1;
            else
                break;
            end
        end
        % Check if the optimized solution is not worse than the initial
        if round(Hstar(end)*(dec/100))/(dec/100) < round(H0*(dec/100))/(dec/100)
            iiih = iiih +1;
        end
%         exitflagh(ii) = exitflag1;
    end
    Heff(ii) = Hstar(end); % Hover efficiency in direction d
    alphastarH(:,ii) = alphastar(:, end); % rotors orientations needed to hover in this mode 
    wstarH(:,ii) = wstar(:,end);% rotors speed needed to hover in this mode 
    
    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(Hstar(end)*dec)/dec <= round(H0*dec)/dec
        worthH = worthH+1;
    end
end
Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Heff = 100*Heff;
Heff = round(Heff*dec/1000)/(dec/1000);
Fnorm = vecnorm(F);
Fmax = max(Fnorm);
Fmin = min(Fnorm);
Mnorm = vecnorm(M);
Mmax = max(Mnorm);
Mmin = min(Mnorm);
Hmax = max(Heff);
Hmin = min(Heff);

worthF = length_D-worthF;
worthM = length_D-worthM;
worthH = length_D-worthH;
% Surface Reconstruction from scattered points cloud
TRI = MyCrustOpen(D_unit.');
[row, ~] = size(TRI);
F_surf =0;
F_vol =0;
M_surf =0;
M_vol =0;

for i = 1:row
    %% Find force space surface and volume:
    AB = F(:,TRI(i, 2))-F(:,TRI(i, 1)); % AB = OB-OA
    AC = F(:,TRI(i, 3))-F(:,TRI(i, 1)); % AC = OC-OA
    AD = [0;0;0]-F(:,TRI(i, 1)); % AD = OD-OA
    F_surf = F_surf + norm(cross(AB,AC))/2;
    F_vol =F_vol + abs(det([AB, AC, AD]))/6;
    %% Find torque space surface and volume:
    AB = M(:,TRI(i, 2))-M(:,TRI(i, 1)); % AB = OB-OA
    AC = M(:,TRI(i, 3))-M(:,TRI(i, 1)); % AC = OC-OA
    AD = [0;0;0]-M(:,TRI(i, 1)); % AD = OD-OA
    M_surf = M_surf + norm(cross(AB,AC))/2;
    M_vol =M_vol + abs(det([AB, AC, AD]))/6;
    %% As a test find the unit sphere surface and volume:
%     AB = D_unit(:,TRI(i, 2))-D_unit(:,TRI(i, 1)); % AB = OB-OA
%     AC = D_unit(:,TRI(i, 3))-D_unit(:,TRI(i, 1)); % AC = OC-OA
%     AD = [0;0;0]-D_unit(:,TRI(i, 1)); % AD = OD-OA
%     D_surf = D_surf + norm(cross(AB,AC))/2;
%     D_vol =D_vol + abs(det([AB, AC, AD]))/6;
end
end

