function [D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, opt_iterations, alphadotmax)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmax, g, step, optim, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_COMPUTE_METRICS computes a lot of metrics for a given design of Quadcopter
%   Design defined by the arms angle (beta & theta) and other parameters

%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
m = Mb+4*Mp;% drone mass [kg]
Ndecimals = 5;
dec = 10.^Ndecimals;
%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb0 = rotz(rad2deg(roll0))*roty(rad2deg(pitch0))*rotz(rad2deg(yaw0)); % Rotation Matrix mapping body frame to inertial frame
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
A_F_static = [sin(beta(1))*cos(theta(1)), sin(theta(1)), -sin(beta(2))*sin(theta(2)), ...
              cos(theta(2)), -sin(beta(3))*cos(theta(3)), -sin(theta(3)), ...
              sin(beta(4))*sin(theta(4)), -cos(theta(4)); ...
              sin(beta(1))*sin(theta(1)), -cos(theta(1)), sin(beta(2))*cos(theta(2)), ...
              sin(theta(2)), -sin(beta(3))*sin(theta(3)), cos(theta(3)), ...
              -sin(beta(4))*cos(theta(4)), -sin(theta(4)); ...
              cos(beta(1)), 0, cos(beta(2)), 0,  cos(beta(3)), 0, cos(beta(4)), 0];

% This static matrix links Fdec to the torque applied by the propellers to
% the drone body 
% M = Ib*wb' = A_M_static*Fdec
A_M_static =[L*sin(theta(1))-km*sin(beta(1))*cos(theta(1))/kf, -L*sin(beta(1))*cos(theta(1))-km*sin(theta(1))/kf, ...
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
F = []; % Matrix containing the maximum force appliable by the design in every direction of D
Feff = []; % Vector containing the efficiency of the maximum forces in every direction of D
M = []; % Matrix containing the maximum torques appliable by the design in every direction of D
Meff = [];% Vector containing the efficiency of the maximum torques in every direction of D
Heff = [];% Vector containing the efficiency of hover mode in every direction of D
C = [];% Vector containing the time to chnge from hover in "normal" pose to hover in every direction of D
number_of_directions = 0; % Counter for the number of direction contained in D
worthF = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal force
worthM = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal torque
worthH = 0; % Counter to quantify the efficiency of the fmincom optimization on hover efficiency
worthC = 0; % Counter to quantify the efficiency of the fmincom optimization on changeability

%% Loop to create the direction matrix D with the desired number of direction:
for i = -1:step:1
    for j = -1:step:1
        for k = -1:step:1% number of directions depends only on step size
            d = [i j k].';
            D = [D d];
        end
    end
end
i0 = find(~vecnorm(D)); 
D(:,i0) = [];% Eliminate the [0; 0; 0] direction
D_unit = D./vecnorm(D); % Create a normalized matrix of direction.
D_unit = round(D_unit*dec)/dec;
[D_unit,ia,ic] = unique(D_unit.', 'stable', 'rows');
D_unit = D_unit.'; % Eliminate redundant directions in normalized D
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

%% Loop to perform the optimizations of the max force, torque, hover efficiency in every direction
for d = D_unit
    number_of_directions = number_of_directions+1;
    
    %% find max force in direction d
    % find initial alpha and n for the optimisation to find the max thrust in direction d
    Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers
    Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
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
    [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta,n0, L, g, Mb, Mp, R, false);
    F0 = m*pdotdot;
    FN0 = norm(F0);
    nstar = n0;
    alphastar = alpha0;
    Fstar = F0;
    FNstar = FN0;
    exitflag1 = 0;
    if optim
        for i = 2:opt_iterations
            % Perform the optimization and find the max thrust in direction d
            [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_thrust(kf, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar(i, :), beta, theta, nstar(:,i), L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            Fstar(:,i) = m*pdotdot;
            FNstar(i) = norm(Fstar(:,i));
            % verify that the found solution satisfies the constraints
            
            if exitflag == -2
                if ~isequal(round(cross(Fstar(:,i),d)*10^2)/10^2,[0 0 0].')
                    nstar(:,i) = nstar(:,i-1);
                    Fstar(:,i) = Fstar(:,i-1);
                    FNstar(i) = FNstar(i-1);
                    break;
                end
                for j =1:4
                    if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                        nstar(:,i) = nstar(:,i-1);
                        Fstar(:,i) = Fstar(:,i-1);
                        FNstar(i) = FNstar(i-1);
                        break;
                    end
                end
                    break;
            end
            % verify that the previous solution is not better
            if round(FNstar(i)*10^4)/10^4 < round(FNstar(i-1)*10^4)/10^4 && i > 2
                    nstar(:,i) = nstar(:,i-1);
                    Fstar(:,i) = Fstar(:,i-1);
                    FNstar(i) = FNstar(i-1);
                    break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(FNstar(i)*10^4)/10^4 == round(FNstar(i-1)*10^4)/10^4 && i > 2
                    break;
            end
        end
        if round(FNstar(end)*10^3)/10^3 < round(FN0*10^3)/10^3
            iiif = iiif +1;
        end
        exitflagf = [exitflagf exitflag1];
    end
    F = [F Fstar(:,end)];% Force produced by the MAV
    Feff = [Feff kf*(nstar(1,end)^2+nstar(2,end)^2+nstar(3,end)^2+nstar(4,end)^2)];
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
    [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta, n0, L, g, Mb, Mp, R, false);
    wbdot = round(dec*wbdot)/dec;
    M0 = Ib*wbdot;
    MN0 = norm(M0);
    nstar = n0;
    alphastar = alpha0;
    Mstar = M0;
    MNstar = MN0;
    exitflag1 = 0;
    if optim
        for i = 2:opt_iterations
            % Perform the optimization and find the max torque in direction d
            [alphastarloop, nstarloop, exitflag] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), d, beta, theta, Display, Algorithm, maxIter,StepTolerance, ConstraintTolerance);
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar(i, :), beta, theta, nstar(:,i), L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar(:,i) = Ib*wbdot;
            MNstar(i) = norm(Mstar(:,i));
            % verify that the found solution satisfies the constraints
            if exitflag == -2
                if ~isequal(round(cross(Mstar(:,i),d)*10^2)/10^2,[0 0 0].')
                    nstar(:,i) = nstar(:,i-1);
                    Mstar(:,i) = Mstar(:,i-1);
                    MNstar(i) = MNstar(i-1); 
                    break;
                end
                for j =1:4
                    if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                        nstar(:,i) = nstar(:,i-1);
                        Mstar(:,i) = Mstar(:,i-1);
                        MNstar(i) = MNstar(i-1); 
                        break;
                    end
                end
                break;
            end
            % verify that the previous solution is not better
            if round(MNstar(i)*10^4)/10^4 < round(MNstar(i-1)*10^4)/10^4 && i > 2
                    nstar(:,i) = nstar(:,i-1);
                    Mstar(:,i) = Mstar(:,i-1);
                    MNstar(i) = MNstar(i-1);
                    break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(MNstar(i)*10^4)/10^4 == round(MNstar(i-1)*10^4)/10^4 && i > 2
                    break;
            end
        end
        if round(MNstar(end)*10^3)/10^3 < round(MN0*10^3)/10^3
            iiim = iiim +1;
        end
        exitflagm = [exitflagm exitflag1];
    end
    M = [M Mstar(:,end)];% Force produced by the MAV
    Meff = [Meff L*kf*(nstar(1,end)^2+nstar(2,end)^2+nstar(3,end)^2+nstar(4,end)^2)];
    if round(MNstar(end)*10^3)/10^3 <= round(MN0*10^3)/10^3
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
    nstar = n0;
    alphastar = alpha0;
    Hstar = H0;
    exitflag1 = 0;
    if optim
        for i = 2:opt_iterations
            % Perform the optimization and find the best hover in direction d
            [alphastarloop, nstarloop, exitflag, TH] = Quadcopter_tilted_arms_min_hover(kf, Fdes, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, false);
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            Hstar(i) = m*g/(kf*norm(nstar(:,i))^2);
            Hstar(i) = round(dec*Hstar(i))/dec;
            % verify that the found solution satisfies the constraints
            if exitflag == -2
                if ~isequal(round(TH*10^2)/10^2,round(Fdes*10^2)/10^2)
                    Hstar(i) = Hstar(i-1); 
                    break;
                end
                for j =1:4
                    if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                        Hstar(i) = Hstar(i-1); 
                        break;
                    end
                end
                break;
            end
            % verify that the previous solution is not better
            if round(Hstar(i)*10^4)/10^4 < round(Hstar(i-1)*10^4)/10^4 && i > 2
                    Hstar(i) = Hstar(i-1);
                    break;
            end
            exitflag1 = exitflag;
            % converged to an optimal solution
            if round(Hstar(i)*10^4)/10^4 == round(Hstar(i-1)*10^4)/10^4 && i > 2
                    break;
            end
        end
        if round(Hstar(end)*10^3)/10^3 < round(H0*10^3)/10^3
            iiih = iiih +1;
        end
        exitflagh = [exitflagh exitflag1];
    end
    Heff = [Heff Hstar(end)]; % Hover efficiency in direction d
    if round(Hstar(end)*10^3)/10^3 <= round(H0*10^3)/10^3
        worthH = worthH+1;
    end
    
    %% find time to change MAV Hover orientation from [0 0 z] to direction d
    % find initial alpha and n for the optimisation to find hover in orientation d
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
    nstar = n0;
    alphastar = alpha0;
    exitflag1 = 0;
    if optim
        for i = 2:opt_iterations
            % Perform the optimization and find the smallest tilting angles allowing hover in direction d
            [alphastarloop, nstarloop, exitflag, TH] = Quadcopter_tilted_arms_min_hover(kf, Fdes, nmin, nmax, alphamin, alphamax, alphastar(i-1,:), nstar(:,i-1), beta, theta, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, true);
            nstar(:,i) = round(dec*nstarloop)/dec;
            alphastar(i, :) = round(dec*alphastarloop)/dec;
            % verify that the found solution satisfies the constraints
            if exitflag == -2
                if ~isequal(round(TH*10^2)/10^2,round(Fdes*10^2)/10^2)
                    nstar(:,i) = nstar(:,i-1);
                    alphastar(i, :) = alphastar(i-1, :);
                    break;
                end
                for j =1:4
                    if nstar(j,i)<nmin || nstar(j,i)>nmax || alphastar(i,j)<alphamin || alphastar(i,j)>alphamax
                        nstar(:,i) = nstar(:,i-1);
                        alphastar(i, :) = alphastar(i-1, :);
                        break;
                    end
                end
                break;
            end
            % verify that the previous solution is not better
            if round(norm(alphastar(i,:))*10^4)/10^4 > round(norm(alphastar(i-1,:))*10^4)/10^4 && i > 2
                    nstar(:,i) = nstar(:,i-1);
                    alphastar(i, :) = alphastar(i-1, :);
                    break;
            end
            exitflag1 =  exitflag;
            % converged to an optimal solution
            if isequal(round(alphastar(i)*10^4)/10^4, round(alphastar(i-1)*10^4)/10^4) && i > 2
                    break;
            end
        end
        if round(norm(alphastar(i,:))*10^3)/10^3 > round(norm(alpha0)*10^3)/10^3
            iiic = iiic +1;
        end
        exitflagc = [exitflagc exitflag1];
    end
    alphastarmax = max(alphastar(end, :));
    tiltingtime = alphastarmax/alphadotmax;
    C = [C, tiltingtime];
    if norm(round(alphastar(i,:)*10^3)/10^3) >= norm(round(alpha0*10^3)/10^3)
        worthC = worthC+1;
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
worthC = number_of_directions-worthC;
% Surface Reconstruction from scattered points cloud
TRI = MyCrustOpen(D_unit.');
[row, column] = size(TRI);
F_surf =0;
F_vol =0;
M_surf =0;
M_vol =0;
% D_surf =0;
% D_vol = 0;
iiif = iiif
iiim = iiim
iiih = iiih
iiic = iiic
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
