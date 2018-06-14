function [wRb, D_unit2, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, length_D, TRI, F_surf, F_vol, M_surf, M_vol] = Mav_compute_metrics_parallel(dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, step, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations)
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
Feff = []; % Vector containing the efficiency of the drone when applying the maximum force in every direction of D
M = []; % Matrix containing the maximum torques appliable by the design in every direction of D
Meff = [];% Vector containing the efficiency of the drone when applying the maximum torque in every direction of D
Heff = []; % Vector containing the efficiency of the drone when hovering with the weight oriented in direction -D
D_unit2 = [];

worthF = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal force
worthM = 0; % Counter to quantify the efficiency of the fmincom optimization on maximal torque
worthH = 0; % Counter to quantify the efficiency of the fmincom optimization on hover efficiency

%% Loop to perform the optimizations of the max force, torque, hover efficiency in every direction
parfor jj = 1:8
    [D_unit1, Heff1, F1, Feff1, M1, Meff1, worthF1, worthM1, worthH1] = Mav_opt_points(jj, wRb, D_unit, length_D, A_F_staticinv, A_M_staticinv, dec, n, beta ,theta, L, kf, km, wmin, wmax, alphamin, alphamax, g, optim, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance, max_iterations);
    Heff = [Heff, Heff1];
    F = [F, F1];
    Feff = [Feff, Feff1];
    M = [M, M1];
    Meff = [Meff, Meff1];
    D_unit2 = [D_unit2, D_unit1];
    worthF = worthF + worthF1;
    worthH = worthH + worthH1;
    worthM = worthM + worthM1;
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

