function [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics(dec, n, beta ,theta, L, Mb, Mp, m, R, kf, km, wmin, wmax, alphamin, alphamax, g, step, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations)
%MAV_COMPUTE_METRICS computes a lot of metrics for a given design of Mav
%   Design defined by the number of arms and their angles (beta & theta) and other parameters

%%%%%%%%%%%% MAV with tilting rotor and tilted arms design optimization%%%%%%%%%%%%

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
% D_unit2 = [];
D = D(:,ia.'); % Eliminate redundant directions in D

%% initialization of the parameters:
[~, length_D] = size(D);

F = zeros(3, length_D); % Matrix containing the maximum force appliable by the design in every direction of D
Feff = zeros(1, length_D); % Vector containing the efficiency of the drone when applying the maximum force in every direction of D
alphastarF = zeros(n, length_D);% Vector containing the optimal tilting angles to obtain the max force in every direction of D
wstarF = zeros(n, length_D); % Vector containing the optimal propeller speed to obtain the max force in every direction of D

M = zeros(3, length_D); % Matrix containing the maximum torques appliable by the design in every direction of D
Meff = zeros(1, length_D);% Vector containing the efficiency of the drone when applying the maximum torque in every direction of D
alphastarM = zeros(n, length_D);% Vector containing the optimal tilting angles to obtain the max torque in every direction of D
wstarM = zeros(n, length_D); % Vector containing the optimal propeller speed to obtain the max torque in every direction of D

Heff = zeros(1, length_D); % Vector containing the efficiency of the drone when hovering with the weight oriented in direction -D 
alphastarH = zeros(n, length_D);% Vector containing the optimal tilting angles to hover in every direction of D
wstarH = zeros(n, length_D); % Vector containing the optimal propeller speed to obtain the max force in every direction of D

number_of_directions = 0; % Counter for the number of direction contained in D

worthF = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal force
worthM = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal torque
worthH = 0; % Counter to quantify the efficiency of the fmincom optimization on hover efficiency

%% Test parameters:
iiif = 0; % incremented if the optimization returns a worse solution than the initial to the maximum force problem in one direction
iiim = 0; % incremented if the optimization returns a worse solution than the initial to the maximum torque problem in one direction
iiih = 0; % incremented if the optimization returns a worse solution than the initial to the optimal hover problem in one direction

exitflagf = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the maximal force in direction d
exitflagm = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the maximal torque in direction d
exitflagh = zeros(1, length_D); % Vector containing all the return flags of fmincom when searching the optimal hover in direction d

% successh = 0;
% successf = 0;
% successm = 0;
% failh = 0;
% failf = 0;
% failm = 0;

%% Loop to perform the optimizations of the max force, torque, hover efficiency in every direction
for d = D_unit
%     D_unit2 = [D_unit2 d];
    number_of_directions = number_of_directions+1; % Count the number of directions
    
    %% First, find the max force in direction d using static matrix
    % Set the initial desired force to be the force to hover in direction d
    Fdes = m*g*d;
    % Loop to find the maximal force appliable by the drone in direction d
    for i = 1:max_iterations
        Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
        
        % Retrieve rotors speeds and orientations from Fdec
        [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
        
        alpha_bound = [alpha0(alpha0>alphamax), alpha0(alpha0<alphamin)];
        w_bound = [w0(w0>wmax), w0(w0<wmin)];
        if isempty(w_bound) && isempty(alpha_bound)
            % Slowly increase Fdes until the obtained alpha0 and w0 does
            % not respect their bounds anymore.
            Fdes = Fdes + 2*d*(n*wmax^2*kf-m*g)/max_iterations; 
        else
            % If alpha0 and w0 does not respect their bounds anymore.
            % Return to the previous Fdes
            Fdes = Fdes - 2*d*(n*wmax^2*kf-m*g)/max_iterations;
            break;
        end
    end
    Fdec = A_F_staticinv*Fdes; % Fdec = inv(Astatic)*Fdes
    
    % Retrieve rotors speeds and orientations from Fdec
    [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
    
    % calculate linear acceleration with this alphastar and nstar
    [m, ~,pdotdot, ~] = Mav_dynamic(n, kf, km, wRb, alpha0, beta, theta,w0, L, g, Mb, Mp, R, dec, false);
    F0 = m*pdotdot;
    FN0 = norm(F0);
    wstar = w0;
    alphastar = alpha0;
    Fstar = F0;
    FNstar = FN0;
    exitflag1 = 0; 
    %% find the max force in direction d using fmincom and static matrix solution as initial solution
    if optim % performs the optimisation only if optim is true
        i = 2;
        while( i < max_iterations) % loop that performs the optimization until the solution is the best possible.
            % Perform the optimization and find the max thrust in direction d
            [alphastarloop, nstarloop, exitflag] = Mav_maximum_force(n, kf, wmin, wmax, alphamin, alphamax, alphastar(:,i-1), wstar(:,i-1), d, beta, theta, wRb, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            
            wstar(:,i) = round(dec*nstarloop)/dec; % rotors speeds after optimization
            alphastar(:, i) = round(dec*alphastarloop)/dec; % rotors orientations after optimization

            % Calculate angular and linear acceleration with this alphastar and nstar
            [m, ~, pdotdot, ~] = Mav_dynamic(n, kf, km, wRb, alphastar(:, i), beta, theta, wstar(:,i), L, g, Mb, Mp, R, dec, false);
            pdotdot = round(dec*pdotdot)/dec;
            Fstar(:,i) = m*pdotdot; % Force applied to the body with the propellers in this 
            FNstar(i) = norm(Fstar(:,i));
            
            % If this solution breaks the constraint: Fstar parallel to d
            % Return to the previous solution and quit the loop
            if ~isequal(round((Fstar(:,i)/norm(Fstar(:,i)))*10^2)/10^2,round(d*10^2)/10^2) 
                alphastar(:,i) = alphastar(:,i-1);
                wstar(:,i) = wstar(:,i-1);
                Fstar(:,i) = Fstar(:,i-1);
                FNstar(i) = FNstar(i-1);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j,i)<alphamin || alphastar(j,i)>alphamax
                    alphastar(:,i) = alphastar(:, i-1);
                    wstar(:,i) = wstar(:,i-1);
                    Fstar(:,i) = Fstar(:,i-1);
                    FNstar(i) = FNstar(i-1);
                    break;
                end
            end
            
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(FNstar(i)*dec)/dec < round(FNstar(i-1)*dec)/dec && i > 3
                alphastar(:, i) = alphastar(:, i-1);
                wstar(:,i) = wstar(:,i-1);
                Fstar(:,i) = Fstar(:,i-1);
                FNstar(i) = FNstar(i-1);
                break;
            end
            exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(FNstar(i)*dec)/dec == round(FNstar(i-2)*dec)/dec
                break;
            end
            i = i+1;
        end
        if round(FNstar(end)*dec)/dec < round(FN0*dec)/dec
            iiif = iiif +1;
        end
        exitflagf(number_of_directions) = exitflag1;
    end
%     [~, lengthF] = size(alphastar);
%     if lengthF > 198
%         kkkkk = lengthF;
%     end
    F(:,number_of_directions) = Fstar(:,end);% Force produced by the MAV
    Feff(number_of_directions) = kf*(norm(wstar(:,end))^2); % Efficiency of this force 
    alphastarF(:,number_of_directions) = alphastar(:,end);% rotors orientations needed to obtain that force
    wstarF(:,number_of_directions) = wstar(:,end);% rotors speed needed to obtain that force
    
    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(FNstar(end)*dec)/dec <= round(FN0*dec)/dec
        worthF = worthF+1;
    end
    
    %% find max torque in direction d using static matrix
    % Set the initial desired torque to be the torque produced if the force
    % to hover was applied at the end of one of the arms
    Mdes = d*(m*g*L);
    % Loop to find the maximal torque appliable by the drone in direction d
    for i = 1:max_iterations
        Fdec = A_M_staticinv*Mdes; % Fdec = inv(Astatic)*Fdes
        
        % Retrieve rotors speeds and orientations from Fdec
        [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);
        
        alpha_bound = [alpha0(alpha0>alphamax), alpha0(alpha0<alphamin)];
        w_bound = [w0(w0>wmax), w0(w0<wmin)];
        if isempty(w_bound) && isempty(alpha_bound)
            % Slowly increase Fdes until the obtained alpha0 and w0 does
            % not respect their bounds anymore.
            Mdes = Mdes + 2*d*(n*L*wmax^2*kf-m*g*L)/max_iterations; 
        else
            % If alpha0 and w0 does not respect their bounds anymore.
            % Return to the previous Fdes
            Mdes = Mdes - 2*d*(n*L*wmax^2*kf-m*g*L)/max_iterations;
            break;
        end
    end
    Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
    
    % Retrieve rotors speeds and orientations from Fdec
    [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec);

    % calculate angular acceleration with this alphastar and nstar
    [~, Ib, ~, wbdot] = Mav_dynamic(n, kf, km, wRb, alpha0, beta, theta, w0, L, g, Mb, Mp, R, dec, false);
    wbdot = round(dec*wbdot)/dec;
    M0 = Ib*wbdot;
    MN0 = norm(M0);
    wstar = w0;
    alphastar = alpha0;
    Mstar = M0;
    MNstar = MN0;
    exitflag1 = 0;
    %% find the max torque in direction d using fmincom and static matrix solution as initial solution
    if optim
        i = 2;
        while( i < opt_iterations)
            % Perform the optimization and find the max torque in direction d
            [alphastarloop, nstarloop, exitflag] = Mav_maximum_torque(n, kf, km, L, wmin, wmax, alphamin, alphamax, alphastar(:,i-1), wstar(:,i-1), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
           
            wstar(:,i) = round(dec*nstarloop)/dec; % optimal rotors speeds
            alphastar(:, i) = round(dec*alphastarloop)/dec; % optimal rotors orientations
            
            % Calculate angular acceleration with this alphastar and nstar
            [~, Ib, ~, wbdot] = Mav_dynamic(n, kf, km, wRb, alphastar(:,i), beta, theta, wstar(:,i), L, g, Mb, Mp, R, dec, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar(:,i) = Ib*wbdot;
            MNstar(i) = norm(Mstar(:,i));
            
            % If this solution breaks the constraint: Mstar parallel to d
            % Return to the previous solution and quit the loop
            if ~isequal(round((Mstar(:,i)/norm(Mstar(:,i)))*10^2)/10^2,round(d*10^2)/10^2)
                wstar(:,i) = wstar(:,i-1);
                alphastar(:, i) = alphastar(:, i-1);
                Mstar(:,i) = Mstar(:,i-1);
                MNstar(i) = MNstar(i-1);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j, i)<alphamin || alphastar(j, i)>alphamax
                    wstar(:,i) = wstar(:,i-1);
                    alphastar(:, i) = alphastar(:, i-1);
                    Mstar(:,i) = Mstar(:,i-1);
                    MNstar(i) = MNstar(i-1);
                    break;
                end
            end
            
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(MNstar(i)*10^2)/10^2 < round(MNstar(i-1)*10^2)/10^2 && i > 3
                    wstar(:,i) = wstar(:,i-1);
                    alphastar(:, i) = alphastar(:, i-1);
                    Mstar(:,i) = Mstar(:,i-1);
                    MNstar(i) = MNstar(i-1);
                    break;
            end
            exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(MNstar(i)*10^4)/10^4 == round(MNstar(i-2)*10^4)/10^4
                break;
            end
            i = i+1;
        end
        
        % Check if the optimized solution is not worse than the initial
        if round(MNstar(end)*10^3)/10^3 < round(MN0*10^3)/10^3
            iiim = iiim +1;
        end
        if ~isempty(m1)
            m3 = Mstar(:,end);
            if round(norm(m1)*10^2)/10^2  > round(norm(m3)*10^2)/10^2
                failm = failm + 1;
            else
                successm = successm + 1;
            end
        end
        exitflagm(number_of_directions) = exitflag1;
    end
%     [lengthM, ~] = size(alphastar);
%     if lengthM > 198
%         lengthM = lengthM;
%     end
    M(:,number_of_directions) = Mstar(:,end);% Torque produced by the MAV
    Meff(number_of_directions) = L*kf*(norm(wstar(:,end))^2); % Efficiency of this torque 
    alphastarM(:,number_of_directions) = alphastar(:, end); % rotors orientations needed to obtain that torque
    wstarM(:,number_of_directions) = wstar(:,end);% rotors speed needed to obtain that torque
    
    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(MNstar(end)*10^2)/10^2 <= round(MN0*10^2)/10^2
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
    wstar = w0;
    alphastar = alpha0;
    Hstar = H0;
    exitflag1 = 0;
    %% find the best hover efficiency in direction d using fmincom and static matrix solution as initial solution
    if optim
        i = 2;
        while( i < opt_iterations)
            
            % Perform the optimization and find the best hover in direction d
            [alphastarloop, nstarloop, exitflag, TH] = Mav_optimize_hover(n, kf, Fdes, wmin, wmax, alphamin, alphamax, alphastar(:, i-1), wstar(:,i-1), beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance);
            wstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(:, i) = round(dec*alphastarloop)/dec;
            Hstar(i) = m*g/(kf*norm(wstar(:,i))^2);
            Hstar(i) = round(dec*Hstar(i))/dec;
            
            % If this solution breaks the constraint: TH = Fdes
            % Return to the previous solution and quit the loop
            if ~isequal(round(TH*10^2)/10^2,round(Fdes*10^2)/10^2)
                alphastar(:, i) = alphastar(:, i-1);
                wstar(:,i) = wstar(:,i-1);
                Hstar(i) = Hstar(i-1);
                break;
            end
            
            % If this solution breaks the constraint: lb < alphastar < ub, lb < wstar <ub
            % Return to the previous solution and quit the loop
            for j =1:n
                if wstar(j,i)<wmin || wstar(j,i)>wmax || alphastar(j, i)<alphamin || alphastar(j, i)>alphamax
                    alphastar(:, i) = alphastar(:, i-1);
                    wstar(:,i) = wstar(:,i-1);
                    Hstar(i) = Hstar(i-1);
                    break;
                end
            end            
           
            % Verify that the previous solution is not better
            % If yes, return to the previous solution and quit the loop
            if round(Hstar(i)*10^3)/10^3 < round(Hstar(i-2)*10^3)/10^3 && i > 5
                    alphastar(:, i) = alphastar(:, i-1);
                    wstar(:,i) = wstar(:,i-1);
                    Hstar(i) = Hstar(i-1);
                    break;
            end
            exitflag1 = exitflag;
            
            % If the loop converged to an optimal solution: quit the loop
            if round(Hstar(i)*10^4)/10^4 == round(Hstar(i-2)*10^4)/10^4
                break;
            end
            i = i+1;
        end
        
        % Check if the optimized solution is not worse than the initial
        if round(Hstar(end)*10^3)/10^3 < round(H0*10^3)/10^3
            iiih = iiih +1;
        end
        exitflagh(number_of_directions) = exitflag1;
    end
%     [lengthH, ~] = size(alphastar);
%     if lengthH > 198
%         lengthH = lengthH;
%     end
    Heff(number_of_directions) = Hstar(end); % Hover efficiency in direction d
    alphastarH(:,number_of_directions) = alphastar(:, end); % rotors orientations needed to hover in this mode 
    wstarH(:,number_of_directions) = wstar(:,end);% rotors speed needed to hover in this mode 
    
    % Test if the solution of the optimization is better than the static
    % matrix one 
    if round(Hstar(end)*10^3)/10^3 <= round(H0*10^3)/10^3
        worthH = worthH+1;
    end
end

Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Heff = 100*Heff;
Heff = round(Heff*dec/100)/(dec/100);
Fnorm = vecnorm(F);
Fmax = max(Fnorm);
Fmin = min(Fnorm);
Mnorm = vecnorm(M);
Mmax = max(Mnorm);
Mmin = min(Mnorm);
Hmax = max(Heff);
Hmin = min(Heff);

worthF = number_of_directions-worthF;
worthM = number_of_directions-worthM;
worthH = number_of_directions-worthH;
% Surface Reconstruction from scattered points cloud
TRI = MyCrustOpen(D_unit.');
[row, ~] = size(TRI);
F_surf =0;
F_vol =0;
M_surf =0;
M_vol =0;

% successh = successh
% successf = successf
% successm = successm
% successc = successc
% 
% failh = failh
% failf = failf
% failm = failm
% failc = failc

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

