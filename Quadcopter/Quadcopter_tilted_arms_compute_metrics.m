function [wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, opt_iterations, alphadotmax)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmax, g, step, optim, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_COMPUTE_METRICS computes a lot of metrics for a given deSign of Quadcopter
%   DeSign defined by the arms angle (beta & theta) and other parameters

%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms deSign optimization%%%%%%%%%%%%
%% Parameters
Ndecimals = 5;
dec = 10.^Ndecimals;
%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb = rotz(rad2deg(yaw0))*roty(rad2deg(pitch0))*rotz(rad2deg(roll0)); % Rotation Matrix mapping body frame to inertial frame
% [m, Ib,pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, [0 0 0 0], beta, theta,[0 0 0 0].', L, g, Mb, Mp, R, false);
%% Re-orientation of the drone so the max force and torque is in direction [0; 0; z]
% if beta(1) ==-beta(3) && beta(1)~=0
%     pitch = beta(3);
% elseif abs(beta(1)) < abs(beta(3))
%     pitch = Sign(beta(3))*(abs(beta(3))-abs(beta(1)))/2;
% elseif abs(beta(3)) < abs(beta(1))
%     pitch = -Sign(beta(1))*(abs(beta(1))-abs(beta(3)))/2;
% else
%     pitch = 0;
% end
% if beta(2) ==-beta(4)&& beta(2)~=0
%     pitch = pitch-pi/2;
%     roll = beta(2)*cos(beta(1));
%     yaw = 0;
%     wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
%     wRb = roty(rad2deg(pi/2))*wRb;
% elseif abs(beta(2)) < abs(beta(4))
%     pitch = pitch + pi/2;
%     roll = cos(beta(1))*Sign(beta(4))*(abs(beta(4))-abs(beta(2)))/2;
%     yaw = 0;
%     wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
%     wRb = roty(rad2deg(-pi/2))*wRb;
% elseif abs(beta(4)) < abs(beta(2))
%     pitch = pitch + pi/2;
%     roll = cos(beta(1))*Sign(beta(2))*(abs(beta(4))-abs(beta(2)))/2;
%     yaw = 0;
%     wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
%     wRb = roty(rad2deg(-pi/2))*wRb;
% else
%     roll = 0;
%     yaw = 0;
%     wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
% end

%% Static matrix
% The static matrix are static allocation matrix that do not depend on the
% rotor orientation and speed.
% (static matrix found using the file:Quadcopter_tilted_arms_Find_Static_Matrix.m)

% Vector containing all the decomposed vertical and horizontal forces:
% Fdec = [kf*cos(alpha(1))*n(1)^2; kf*sin(alpha(1))*n(1)^2; kf*cos(alpha(2))*n(2)^2; 
%         kf*sin(alpha(2))*n(2)^2; kf*cos(alpha(3))*n(3)^2; kf*sin(alpha(3))*n(3)^2; 
%         kf*cos(alpha(4))*n(4)^2; kf*sin(alpha(4))*n(4)^2];

% This static matrix links Fdec to the force applied by the propellers to
% the drone body 
% F = m*p'' = A_F_static*Fdec
A_F_static = wRb*[sin(beta(1))*cos(theta(1)), sin(theta(1)), -sin(beta(2))*sin(theta(2)), ...
              cos(theta(2)), -sin(beta(3))*cos(theta(3)), -sin(theta(3)), ...
              sin(beta(4))*sin(theta(4)), -cos(theta(4)); ...
              sin(beta(1))*sin(theta(1)), -cos(theta(1)), sin(beta(2))*cos(theta(2)), ...
              sin(theta(2)), -sin(beta(3))*sin(theta(3)), cos(theta(3)), ...
              -sin(beta(4))*cos(theta(4)), -sin(theta(4)); ...
              cos(beta(1)), 0, cos(beta(2)), 0,  cos(beta(3)), 0, cos(beta(4)), 0];

% This static matrix links Fdec to the torque applied by the propellers to
% the drone body 
% M = Ib*wb' = A_M_static*Fdec
A_M_static = [L*sin(theta(1))-km*sin(beta(1))*cos(theta(1))/kf, -L*sin(beta(1))*cos(theta(1))-km*sin(theta(1))/kf, ...
             L*cos(theta(2))-km*sin(beta(2))*sin(theta(2))/kf, L*sin(beta(2))*sin(theta(2))+km*cos(theta(2))/kf, ...
             -L*sin(theta(3))+km*sin(beta(3))*cos(theta(3))/kf, L*sin(beta(3))*cos(theta(3))+km*sin(theta(3))/kf, ...
             -L*cos(theta(4))-km*sin(beta(4))*sin(theta(4))/kf, -L*sin(beta(4))*sin(theta(4))-km*cos(theta(4))/kf; ...
             -L*cos(theta(1))-km*sin(beta(1))*sin(theta(1))/kf, -L*sin(beta(1))*sin(theta(1))+km*cos(theta(1))/kf, ...
             L*sin(theta(2))+km*sin(beta(2))*cos(theta(2))/kf, -L*sin(beta(2))*cos(theta(2))+km*sin(theta(2))/kf, ...
             L*cos(theta(3))+km*sin(beta(3))*sin(theta(3))/kf, L*sin(beta(3))*sin(theta(3))-km*cos(theta(3))/kf, ...
             -L*sin(theta(4))-km*sin(beta(4))*cos(theta(4))/kf, L*sin(beta(4))*cos(theta(4))-km*sin(theta(4))/kf; ...
             -km*cos(beta(1))/kf, -L*cos(beta(1)), km*cos(beta(2))/kf, -L*cos(beta(2)), -km*cos(beta(3))/kf, ...
             -L*cos(beta(3)), km*cos(beta(4))/kf, -L*cos(beta(4))];

% The Moore-Penrose pseudo inverse of the static matrices allow to find
% Fdec from a desired force or torque applied on the drone.
% The rotor orientation and speed can then be deduced from Fdec

% => Fdec = inv(A_F_static)*Fdes
A_F_staticinv = pinv(A_F_static);

% => Fdec = inv(A_M_static)*Mdes
A_M_staticinv = pinv(A_M_static);

%% initialization of the optimization:
D = []; % Matrix containing every direction for which we want to compute the metrics
F = []; % Matrix containing the maximum force appliable by the deSign in every direction of D
Feff = []; % Vector containing the efficiency of the maximum forces in every direction of D
alphastarF = [];
nstarF = []; 
alphastarM = [];
nstarM = []; 
alphastarH = [];
nstarH = []; 
alphastarC = [];
nstarC = []; 
M = []; % Matrix containing the maximum torques appliable by the deSign in every direction of D
Meff = [];% Vector containing the efficiency of the maximum torques in every direction of D
Heff = [];% Vector containing the efficiency of hover mode in every direction of D
C = [];% Vector containing the time to chnge from hover in "normal" pose to hover in every direction of D
number_of_directions = 0; % Counter for the number of direction contained in D
worthF = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal force
worthM = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal torque
worthH = 0; % Counter to quantify the efficiency of the fmincom optimization on hover efficiency
worthC = 0; % Counter to quantify the efficiency of the fmincom optimization on changeability

%% Loop to create the direction matrix D with the desired number of direction:
for i = 1:-step:-1
    for j = 1:-step:-1
        for k = 1:-step:-1% number of directions depends only on step size
            d = [i j k].';
            D = [D d];
        end
    end
end
i0 = ~vecnorm(D); 
D(:,i0) = [];% Eliminate the [0; 0; 0] direction
D_unit = D./vecnorm(D); % Create a normalized matrix of direction.
D_unit = round(D_unit*dec)/dec;
[D_unit,ia,~] = unique(D_unit.', 'stable', 'rows');
D_unit = D_unit.'; % Eliminate redundant directions in normalized D
D_unit2 = [];
D = D(:,ia.');% Eliminate redundant directions in D

%% Tests:
iiih = 0;
iiif = 0;
iiim = 0;
iiic = 0;
exitflagh = [];
exitflagf = [];
exitflagm = [];
exitflagc = [];

successh = 0;
successf = 0;
successm = 0;
successc = 0;

failh = 0;
failf = 0;
failm = 0;
failc = 0;

%% Loop to perform the optimizations of the max force, torque, hover efficiency in every direction
for d = D_unit
    D_unit2 = [D_unit2 d];
    number_of_directions = number_of_directions+1; % Count the number of directions
    
    %% find the max force in direction d using static matrix
    Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers in direction d
    
    Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
    
    % Inverse substitution :
    %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
    n0 = n0/vecnorm(n0);
    n0 = nmax^2*n0/max(n0); % nstar <= nmax
    n0(n0<nmin^2) = nmin^2; % nmin <= nstar
    n0 = sqrt(n0);
    alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
    alpha0 = round(dec*alpha0)/dec; % rotors speeds 
    n0 = round(dec*n0)/dec; % rotors orientations after optimization
    
    % calculate angular and linear acceleration with this alphastar and nstar
    [m, ~,pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha0, beta, theta,n0, L, g, Mb, Mp, R, false);
    F0 = m*pdotdot;
    FN0 = norm(F0);
    nstar = [n0 n0];
    alphastar = [alpha0; alpha0];
    Fstar = [F0 F0];
    FNstar = [FN0 FN0];
    exitflag1 = 0; 
    
    %% find the max force in direction d using fmincom and static matrix solution as initial solution
    if optim % performs the optimisation only if optim is true
        f1 = [];
        i = 3;
        while( i < opt_iterations) % loop that performs the optimization until the solution is the best possible.
            % Perform the optimization and find the max thrust in direction d
            if FNstar(i-1) >0
                [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_thrust(kf, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), d, beta, theta, wRb, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            else
                [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_thrust(kf, nmin, nmax, alphamin, alphamax, alphastar(i-2,:), nstar(:,i-2), d, beta, theta, wRb, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            end
            nstar(:,i) = round(dec*nstarloop)/dec; % rotors speeds after optimization
            alphastar(i, :) = round(dec*alphastarloop)/dec; % rotors orientations after optimization

            % calculate angular and linear acceleration with this alphastar and nstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar(i, :), beta, theta, nstar(:,i), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Fstar(:,i) = m*pdotdot; % Force applied to the body with the propellers in this 
            FNstar(i) = norm(Fstar(:,i));
            % if this solution breaks the constraint Fstar // d
            if ~isequal(round((Fstar(:,i)/norm(Fstar(:,i)))*10^2)/10^2,round(d*10^2)/10^2) 
                alphastar(i,:) = alphastar(i-2,:);
                nstar(:,i) = nstar(:,i-2);
                Fstar(:,i) = Fstar(:,i-2);
                FNstar(i) = FNstar(i-2);
                break;
            end
            for j =1:4
                if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                    alphastar(i,:) = alphastar(i-2,:);
                    nstar(:,i) = nstar(:,i-2);
                    Fstar(:,i) = Fstar(:,i-2);
                    FNstar(i) = FNstar(i-2);
                    break;
                end
            end
            % verify that the previous solution is not better
            if round(FNstar(i)*10^4)/10^4 < round(FNstar(i-2)*10^4)/10^4 && i > 5
                alphastar(i, :) = alphastar(i-2, :);
                nstar(:,i) = nstar(:,i-2);
                Fstar(:,i) = Fstar(:,i-2);
                FNstar(i) = FNstar(i-2);
                break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(FNstar(i)*10^4)/10^4 == round(FNstar(i-2)*10^4)/10^4
                %As the design is symetric, look for better solutions in
                %oppoit directions
                dd = find(D_unit2(1,:) == - d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(2) = -Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(3) = -Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            alpha(4) = -Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % verify that this solution is better than the previous and does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(2) = Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(3) = -Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            alpha(4) = Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(2) = -Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(3) = Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            alpha(4) = -Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha = -alphastarF(dd,:);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -alphastarF(dd,1);
                            alpha(2) = alphastarF(dd,2);
                            alpha(3) = -alphastarF(dd,3);
                            alpha(4) = alphastarF(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(2) = Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(3) = Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            alpha(4) = Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = alphastarF(dd,1);
                            alpha(2) = -alphastarF(dd,2);
                            alpha(3) = alphastarF(dd,3);
                            alpha(4) = -alphastarF(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarF(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarF(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarF(:,dd);
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                %%%Opposit
                dd = find(D_unit2(1,:) == - d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(2) = -Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(3) = -Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            alpha(4) = -Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(2) = Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(3) = -Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            alpha(4) = Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(2) = -Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(3) = Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            alpha(4) = -Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -alphastarF(dd, 2);
                            alpha(2) = -alphastarF(dd, 1);
                            alpha(3) = -alphastarF(dd, 4);
                            alpha(4) = -alphastarF(dd, 3);
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -alphastarF(dd, 2);
                            alpha(2) = alphastarF(dd, 1);
                            alpha(3) = -alphastarF(dd, 4);
                            alpha(4) = alphastarF(dd, 3);
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = Sign(alphastarF(dd, 2))*(pi - abs(alphastarF(dd, 2)));
                            alpha(2) = Sign(alphastarF(dd, 1))*(pi - abs(alphastarF(dd, 1)));
                            alpha(3) = Sign(alphastarF(dd, 4))*(pi - abs(alphastarF(dd, 4)));
                            alpha(4) = Sign(alphastarF(dd, 3))*(pi - abs(alphastarF(dd, 3)));
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = alphastarF(dd, 2);
                            alpha(2) = -alphastarF(dd, 1);
                            alpha(3) = alphastarF(dd, 4);
                            alpha(4) = -alphastarF(dd, 3);
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(FNstar(i)*10^3)/10^3 < round(norm(F(:,dd))*10^3)/10^3
                            f1 = Fstar(:,i);
                            alpha(1) = -alphastarF(dd, 2);
                            alpha(2) = -alphastarF(dd, 1);
                            alpha(3) = -alphastarF(dd, 4);
                            alpha(4) = -alphastarF(dd, 3);
                            n(1) = nstarF(2,dd);
                            n(2) = nstarF(1,dd);
                            n(3) = nstarF(4,dd);
                            n(4) = nstarF(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Ftest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Fstar // d
                            if isequal(round((Ftest/norm(Ftest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Ftest)*10^2)/10^2 > round(FNstar(i)*10^2)/10^2
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Fstar(:,i) = Ftest;
                                    FNstar(i) = norm(Ftest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                break;
            end
            i = i+1;
            if i <= opt_iterations-2
                alphastar(i,:) = [0 0 0 0];
                nstar(:,i) = [0 0 0 0].';
                Fstar(:,i) = [0 0 0].';
                FNstar(i) = 0;
            end
            i = i+1;
        end
        if ~isempty(f1)
            f3 = Fstar(:,end);
            if round(norm(f1)*10^2)/10^2 > round(norm(f3)*10^2)/10^2
                failf = failf + 1;
            else
                successf = successf + 1;
            end
        end
        if round(FNstar(end)*10^3)/10^3 < round(FN0*10^3)/10^3
            iiif = iiif +1;
        end
        exitflagf = [exitflagf exitflag1];
    end
    [lengthF, ~] = size(alphastar);
    if lengthF > 198
        lengthF = lengthF;
    end
    F = [F Fstar(:,end)];% Force produced by the MAV
    Feff = [Feff kf*(nstar(1,end)^2+nstar(2,end)^2+nstar(3,end)^2+nstar(4,end)^2)];
    alphastarF = [alphastarF; alphastar(end,:)];
    nstarF = [nstarF nstar(:,end)];
    if round(FNstar(end)*10^3)/10^3 <= round(FN0*10^3)/10^3
        worthF = worthF+1;
    end
    
    %% find max torque in direction d
    % find initial alpha and n for the optimisation to find the max torque in direction d
    Mdes = d*(4*L*nmax^2*kf); % set desired torque to be equal to the maximal torque of the four propellers
    Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
    % Inverse substitution :
    %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
    n0 = n0/vecnorm(n0);
    n0 = nmax^2*n0/max(n0); % 0 <= nstar <= nmax
    n0(n0<nmin^2) = nmin^2;
    n0 = sqrt(n0);
    alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
    alpha0 = round(dec*alpha0)/dec;
    n0 = round(dec*n0)/dec;
    % calculate angular and linear acceleration with this alphastar and nstar
    [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha0, beta, theta, n0, L, g, Mb, Mp, R, false);
    wbdot = round(dec*wbdot)/dec;
    M0 = Ib*wbdot;
    MN0 = norm(M0);
    nstar = [n0 n0];
    alphastar = [alpha0; alpha0];
    Mstar = [M0 M0];
    MNstar = [MN0 MN0];
    exitflag1 = 0;
    if optim
        m1 =[];
        i = 3;
        while( i < opt_iterations)
            % Perform the optimization and find the max torque in direction d
            if MNstar(i-1) >0
                [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            else
                [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmin, nmax, alphamin, alphamax, alphastar(i-2,:), nstar(:,i-2), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            end            
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar(i, :), beta, theta, nstar(:,i), L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar(:,i) = Ib*wbdot;
            MNstar(i) = norm(Mstar(:,i));
            % verify that the found solution satisfies the constraints
            if ~isequal(round((Mstar(:,i)/norm(Mstar(:,i)))*10^2)/10^2,round(d*10^2)/10^2)
                nstar(:,i) = nstar(:,i-2);
                alphastar(i, :) = alphastar(i-2, :);
                Mstar(:,i) = Mstar(:,i-2);
                MNstar(i) = MNstar(i-2);
                break;
            end
            for j =1:4
                if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                    nstar(:,i) = nstar(:,i-2);
                    alphastar(i, :) = alphastar(i-2, :);
                    Mstar(:,i) = Mstar(:,i-2);
                    MNstar(i) = MNstar(i-2);
                    break;
                end
            end
            % verify that the previous solution is not better
            if round(MNstar(i)*10^2)/10^2 < round(MNstar(i-2)*10^2)/10^2 && i > 5
                    nstar(:,i) = nstar(:,i-2);
                    alphastar(i, :) = alphastar(i-2, :);
                    Mstar(:,i) = Mstar(:,i-2);
                    MNstar(i) = MNstar(i-2);
                    break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(MNstar(i)*10^4)/10^4 == round(MNstar(i-2)*10^4)/10^4
                %As the design is symetric, look for better solutions in
                %oppoit directions
                dd = find(D_unit2(1,:) == - d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(2) = -Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(3) = -Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            alpha(4) = -Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(2) = Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(3) = -Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            alpha(4) = Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(2) = -Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(3) = Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            alpha(4) = -Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(2));
                DD = DD(:,dd2);
                dd3  = find(DD(3,:) == d(3));
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha = -alphastarM(dd,:);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -alphastarM(dd,1);
                            alpha(2) = alphastarM(dd,2);
                            alpha(3) = -alphastarM(dd,3);
                            alpha(4) = alphastarM(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(2) = Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(3) = Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            alpha(4) = Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = alphastarM(dd,1);
                            alpha(2) = -alphastarM(dd,2);
                            alpha(3) = alphastarM(dd,3);
                            alpha(4) = -alphastarM(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarM(:,dd), L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarM(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarM(:,dd);
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                %%%Opposit
                dd = find(D_unit2(1,:) == - d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(2) = -Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(3) = -Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            alpha(4) = -Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(2) = Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(3) = -Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            alpha(4) = Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(2) = -Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(3) = Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            alpha(4) = -Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -alphastarM(dd, 2);
                            alpha(2) = -alphastarM(dd, 1);
                            alpha(3) = -alphastarM(dd, 4);
                            alpha(4) = -alphastarM(dd, 3);
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -alphastarM(dd, 2);
                            alpha(2) = alphastarM(dd, 1);
                            alpha(3) = -alphastarM(dd, 4);
                            alpha(4) = alphastarM(dd, 3);
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = Sign(alphastarM(dd, 2))*(pi - abs(alphastarM(dd, 2)));
                            alpha(2) = Sign(alphastarM(dd, 1))*(pi - abs(alphastarM(dd, 1)));
                            alpha(3) = Sign(alphastarM(dd, 4))*(pi - abs(alphastarM(dd, 4)));
                            alpha(4) = Sign(alphastarM(dd, 3))*(pi - abs(alphastarM(dd, 3)));
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = find(DD(3,:) == d(3));
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = alphastarM(dd, 2);
                            alpha(2) = -alphastarM(dd, 1);
                            alpha(3) = alphastarM(dd, 4);
                            alpha(4) = -alphastarM(dd, 3);
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(MNstar(i)*10^3)/10^3 < round(norm(M(:,dd))*10^3)/10^3
                            m1 = Mstar(:,i);
                            alpha(1) = -alphastarM(dd, 2);
                            alpha(2) = -alphastarM(dd, 1);
                            alpha(3) = -alphastarM(dd, 4);
                            alpha(4) = -alphastarM(dd, 3);
                            n(1) = nstarM(2,dd);
                            n(2) = nstarM(1,dd);
                            n(3) = nstarM(4,dd);
                            n(4) = nstarM(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [~, Ib, ~, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            wbdot = round(dec*wbdot)/dec;
                            Mtest = Ib*wbdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round((Mtest/norm(Mtest))*10^2)/10^2,round(d*10^2)/10^2) && i < opt_iterations-2 && round(norm(Mtest)*10^2)/10^2 > round(MNstar(i)*10^2)/10^2 
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Mstar(:,i) = Mtest;
                                    MNstar(i) = norm(Mtest);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                break;
            end
            i = i+1;
            if i <= opt_iterations-2
                alphastar(i,:) = [0 0 0 0];
                nstar(:,i) = [0 0 0 0].';
                Mstar(:,i) = [0 0 0].';
                MNstar(i) = 0;
            end
            i = i+1;
        end
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
        
        exitflagm = [exitflagm exitflag1];
    end
    [lengthM, ~] = size(alphastar);
    if lengthM > 198
        lengthM = lengthM;
    end
    M = [M Mstar(:,end)];% Force produced by the MAV
    Meff = [Meff L*kf*(nstar(1,end)^2+nstar(2,end)^2+nstar(3,end)^2+nstar(4,end)^2)];
    alphastarM = [alphastarM; alphastar(end,:)];
    nstarM = [nstarM nstar(:,end)];
    if round(MNstar(end)*10^2)/10^2 <= round(MN0*10^2)/10^2
        worthM = worthM+1;
    end

    %% find hover efficiency in direction d
    % find initial alpha and n for the optimisation to find the best hover in orientation d
    Fdes = m*g*d;
    Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
    % Inverse substitution :
    %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
    n0(n0<nmin^2) = nmin^2;
    n0 = sqrt(n0);
    alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
    n0 = round(dec*n0)/dec;
    alpha0 = round(dec*alpha0)/dec;
    H0 = m*g/(kf*norm(n0)^2);
    H0 = round(dec*H0)/dec;
    nstar = [n0 n0];
    alphastar = [alpha0; alpha0];
    Hstar = [H0 H0];
    exitflag1 = 0;
    if optim
        h1 = [];
        i = 3;
        while( i < opt_iterations)
            % Perform the optimization and find the best hover in direction d
            if Hstar(i-1) >0
                [alphastarloop, nstarloop, exitflag, TH] = Quadcopter_tilted_arms_min_hover(kf, Fdes, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, false);
            else
                [alphastarloop, nstarloop, exitflag, TH] = Quadcopter_tilted_arms_min_hover(kf, Fdes, nmin, nmax, alphamin, alphamax, alphastar(i-2,:), nstar(:,i-2), beta, theta, wRb, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, false);
            end
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            Hstar(i) = m*g/(kf*norm(nstar(:,i))^2);
            Hstar(i) = round(dec*Hstar(i))/dec;
            % verify that the found solution satisfies the constraints
            if ~isequal(round(TH*10^2)/10^2,round(Fdes*10^2)/10^2)
                alphastar(i, :) = alphastar(i-2, :);
                nstar(:,i) = nstar(:,i-2);
                Hstar(i) = Hstar(i-2);
                break;
            end
            for j =1:4
                if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                    alphastar(i, :) = alphastar(i-2, :);
                    nstar(:,i) = nstar(:,i-2);
                    Hstar(i) = Hstar(i-2);
                    break;
                end
            end
            % verify that the previous solution is not better
            if round(Hstar(i)*10^3)/10^3 < round(Hstar(i-2)*10^3)/10^3 && i > 5
                    alphastar(i, :) = alphastar(i-2, :);
                    nstar(:,i) = nstar(:,i-2);
                    Hstar(i) = Hstar(i-2);
                    break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(Hstar(i)*10^4)/10^4 == round(Hstar(i-2)*10^4)/10^4
                %As the design is symetric, look for better solutions in
                %oppoit directions
                dd = find(D_unit2(1,:) == - d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(2) = -Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(3) = -Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            alpha(4) = -Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(2) = Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(3) = -Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            alpha(4) = Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(2) = -Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(3) = Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            alpha(4) = -Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha = -alphastarH(dd,:);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -alphastarH(dd,1);
                            alpha(2) = alphastarH(dd,2);
                            alpha(3) = -alphastarH(dd,3);
                            alpha(4) = alphastarH(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(2) = Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(3) = Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            alpha(4) = Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(1));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(2));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = alphastarH(dd,1);
                            alpha(2) = -alphastarH(dd,2);
                            alpha(3) = alphastarH(dd,3);
                            alpha(4) = -alphastarH(dd,4);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, nstarH(:,dd), L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(nstarH(:,dd)*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = nstarH(:,dd);
                                    Hstar(i) = m*g/(kf*norm(nstarH(:,dd))^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                %%%Opposit
                dd = find(D_unit2(1,:) == - d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(2) = -Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(3) = -Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            alpha(4) = -Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(2) = Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(3) = -Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            alpha(4) = Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(2) = -Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(3) = Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            alpha(4) = -Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == -d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -alphastarH(dd, 2);
                            alpha(2) = -alphastarH(dd, 1);
                            alpha(3) = -alphastarH(dd, 4);
                            alpha(4) = -alphastarH(dd, 3);
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == - d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -alphastarH(dd, 2);
                            alpha(2) = alphastarH(dd, 1);
                            alpha(3) = -alphastarH(dd, 4);
                            alpha(4) = alphastarH(dd, 3);
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == - d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = Sign(alphastarH(dd, 2))*(pi - abs(alphastarH(dd, 2)));
                            alpha(2) = Sign(alphastarH(dd, 1))*(pi - abs(alphastarH(dd, 1)));
                            alpha(3) = Sign(alphastarH(dd, 4))*(pi - abs(alphastarH(dd, 4)));
                            alpha(4) = Sign(alphastarH(dd, 3))*(pi - abs(alphastarH(dd, 3)));
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == -d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = alphastarH(dd, 2);
                            alpha(2) = -alphastarH(dd, 1);
                            alpha(3) = alphastarH(dd, 4);
                            alpha(4) = -alphastarH(dd, 3);
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                dd = find(D_unit2(1,:) == d(2));
                DD = D_unit2(:,dd);
                dd2 = find(DD(2,:) == d(1));
                DD = DD(:,dd2);
                dd3  = DD(3,:) == d(3);
                dd = dd(dd2(dd3));
                if ~isequal(D_unit2(:,end), D_unit2(:,dd))
                    if ~isempty(dd)
                        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
                            h1 = Hstar(:,i);
                            alpha(1) = -alphastarH(dd, 2);
                            alpha(2) = -alphastarH(dd, 1);
                            alpha(3) = -alphastarH(dd, 4);
                            alpha(4) = -alphastarH(dd, 3);
                            n(1) = nstarH(2,dd);
                            n(2) = nstarH(1,dd);
                            n(3) = nstarH(4,dd);
                            n(4) = nstarH(3,dd);
                            % calculate angular and linear acceleration with this alphastar and nstar
                            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
                            pdotdot = round(dec*pdotdot)/dec;
                            Htest = m*pdotdot; % Force applied to the body with the propellers in this
                            % if this solution does not break the constraint Mstar // d
                            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                                test = [];
                                [rows, ~] = size(alphastar);
                                for j = 1:rows
                                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(nstar(:, j)*10^3)/10^3)
                                        test = [test true];
                                    end
                                end
                                if isempty(test)
                                    i = i+1;
                                    alphastar(i,:) = alpha;
                                    nstar(:,i) = n;
                                    Hstar(i) = m*g/(kf*norm(n)^2);
                                    i = i+1;
                                    continue;
                                end
                            end
                        end
                    end
                end
                break;
            end
            i = i+1;
            if i <= opt_iterations-2
                alphastar(i,:) = [0 0 0 0];
                nstar(:,i) = [0 0 0 0].';
                Hstar(i) = 0;
            end
            i = i+1;
        end
        if ~isempty(h1)
            h3 = Hstar(:,end);
            if round(norm(h1)*10^2)/10^2 > round(norm(h3)*10^2)/10^2
                failh = failh + 1;
            else
                successh = successh + 1;
            end
        end
        if round(Hstar(end)*10^3)/10^3 < round(H0*10^3)/10^3
            iiih = iiih +1;
        end
        exitflagh = [exitflagh exitflag1];
    end
    [lengthH, ~] = size(alphastar);
    if lengthH > 198
        lengthH = lengthH;
    end
    Heff = [Heff Hstar(end)]; % Hover efficiency in direction d
    alphastarH = [alphastarH; alphastar(end,:)];
    nstarH = [nstarH nstar(:,end)];
    if round(Hstar(end)*10^3)/10^3 <= round(H0*10^3)/10^3
        worthH = worthH+1;
    end
end
% F = wRb.'*F;
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
worthC = number_of_directions-worthC;
% Surface Reconstruction from scattered points cloud
TRI = MyCrustOpen(D_unit.');
[row, ~] = size(TRI);
F_surf =0;
F_vol =0;
M_surf =0;
M_vol =0;
% D_surf =0;
% D_vol = 0;

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
