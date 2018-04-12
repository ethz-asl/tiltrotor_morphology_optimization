function [D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, worthF, worthM, worthH, errF, errM, errH, ii] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmax, g, step, optim, Display, Algorithm, maxIter)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmax, g, step, optim, Display, Algorithm, maxIter)
%QUADCOPTER_TILTED_ARMS_COMPUTE_METRICS compute the maximal thrusts and torques in every direction for an
%                                       arbitrary design of quadcopter 
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
wRb0 = rotz(rad2deg(roll0))*roty(rad2deg(pitch0))*rotz(rad2deg(yaw0));
%% Static matrix definition
% Fdec = [kf*cos(alpha(1))*n(1)^2; kf*sin(alpha(1))*n(1)^2; kf*cos(alpha(2))*n(2)^2; 
%         kf*sin(alpha(2))*n(2)^2; kf*cos(alpha(3))*n(3)^2; kf*sin(alpha(3))*n(3)^2; 
%         kf*cos(alpha(4))*n(4)^2; kf*sin(alpha(4))*n(4)^2];

% F = m*p'' = A_F_static*Fdec
% M = Ib*wb' = A_M_static*Fdec

% These static matrix are found using the file: Quadcopter_tilted_arms_Find_Static_Matrix.m
A_F_static = [sin(beta(1))*cos(theta(1)), sin(theta(1)), -sin(beta(2))*sin(theta(2)), ...
              cos(theta(2)), -sin(beta(3))*cos(theta(3)), -sin(theta(3)), ...
              sin(beta(4))*sin(theta(4)), -cos(theta(4)); ...
              sin(beta(1))*sin(theta(1)), -cos(theta(1)), sin(beta(2))*cos(theta(2)), ...
              sin(theta(2)), -sin(beta(3))*sin(theta(3)), cos(theta(3)), ...
              -sin(beta(4))*cos(theta(4)), -sin(theta(4)); ...
              cos(beta(1)), 0, cos(beta(2)), 0,  cos(beta(3)), 0, cos(beta(4)), 0];


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


% => Fdec = inv(A_F_static)*Fdes
A_F_staticinv = pinv(A_F_static);

% => Fdec = inv(A_M_static)*Mdes
A_M_staticinv = pinv(A_M_static);

%% Optimization of alpha and n 
% initialization:
d = [0 0 0].';
D = d;
F = [0 0 0].';
Feff = 0;
M = zeros(3,1);
Meff = 0;
Heff = 0;

% Test
ii = 0;
errF = 0;
errM = 0;
errH = 0;
worthF = 0;
worthM = 0;
worthH = 0;
% Star = [];
% FH = zeros(3,1);
% StarOpt = zeros(10,1);
% Fcheck= [0 0 0].';
% MM = zeros(3,1);

%% Loop to compute the optimal Force in "any" directions (neglecting gravity):
for i = -1:step:1
    for j = -1:step:1
        for k = -1:step:1
            d = [i j k].';
            D = [D d];
            d0 = [0 0 0].';
            if isequal(d,d0)
                d = [0;0;0];
            else
                d = d./vecnorm(d);
            end
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
            n0(n0<0) = 0; 
            n0 = sqrt(n0);
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alpha0 = round(dec*alpha0)/dec;
            n0 = round(dec*n0)/dec;
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta,n0, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            F0 = m*pdotdot;
            FN0 = norm(F0);
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max thrust in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter);
                    alphastar = round(dec*alphastar)/dec;
                    nstar = round(dec*nstar)/dec;
                    %test
                    if isequal(alphastar, alpha0) && isequal(nstar, n0)
                        worthF = worthF+1;
                    end
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta,nstar, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Fstar = m*pdotdot;
            FNstar = norm(Fstar);
            if FNstar < FN0
                F = [F F0];% Force produced by the MAV
                Feff = [Feff kf*n0.'*n0];
                errF = errF+1;
            else
                F = [F Fstar];% Force produced by the MAV
                Feff = [Feff kf*nstar.'*nstar];
            end

%             % Double check            
%             Fdec_check = A_F_staticinv*(m*pdotdot);
%             nstar1 = [1/kf*sqrt(Fdec_check(1)^2 + Fdec_check(2)^2); 1/kf*sqrt(Fdec_check(3)^2 + Fdec_check(4)^2); 1/kf*sqrt(Fdec_check(5)^2 + Fdec_check(6)^2); 1/kf*sqrt(Fdec_check(7)^2 + Fdec_check(8)^2)];
%             nstar1 = sqrt(nstar1);
%             nstar1 = round(dec*nstar1)/dec;
%             alphastar1 = [atan2(Fdec_check(2),Fdec_check(1)) atan2(Fdec_check(4),Fdec_check(3)) atan2(Fdec_check(6),Fdec_check(5)) atan2(Fdec_check(8),Fdec_check(7))];
%             alphastar1 = round(dec*alphastar1)/dec;
%             if norm(alphastar1) < 4*pi
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar1, beta, theta,nstar1, L, g, Mb, Mp, R, false);
%                 pdotdot = round(dec*pdotdot)/dec;
%             else
%                 pdotdot =[0 0 0].';
%             end
%             Fcheck = [Fcheck m*pdotdot];
%             star = [nstar;alphastar.';nstar1;alphastar1.'];
%             Star = [Star star];
            
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
            n0(n0<0) = 0; 
            n0 = sqrt(n0);
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alpha0 = round(dec*alpha0)/dec;
            n0 = round(dec*n0)/dec;
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta, n0, L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            M0 = Ib*wbdot;
            MN0 = norm(M0);
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max torque in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, alpha0, n0, d, beta, theta, Display, Algorithm, maxIter);
                    alphastar = round(dec*alphastar)/dec;
                    nstar = round(dec*nstar)/dec;
                    % Test
                    if isequal(alphastar, alpha0) && isequal(nstar, n0)
                        worthM = worthM+1;
                    end
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta, nstar, L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar = Ib*wbdot;
            MNstar = norm(Mstar);
            if MNstar < MN0
                M = [M M0];% Force produced by the MAV
                Meff = [Meff L*kf*n0.'*n0];
                errM = errM+1;
            else
                M = [M Mstar];% Force produced by the MAV
                Meff = [Meff L*kf*nstar.'*nstar];
            end
            % Test
%             star = [nstar;alphastar.'];
%             Star = [Star, star];
            %% find hover efficiency in direction d 
            if ~isequal(d, [0 0 0].')
                roll = 0;
                pitch = acos(d(3));
                yaw = atan2(d(2),d(1));
                wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
                % find initial alpha and n for the optimisation to find the max torque in direction d
                Fdes = m*g*d;
                Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
                % Inverse substitution : 
                %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
                %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
                n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
                n0 = sqrt(n0);
                alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
                n0 = round(dec*n0)/dec;
                alpha0 = round(dec*alpha0)/dec;
                %Test
%                 star = [round(n0); round(rad2deg(alpha0)).'; 0; round(100*m*g/(kf*norm(n0)^2))];
%                 Star = [Star star];
                if optim
                    % Perform the optimization and find the best hover in direction d 
                    [alphastar, nstar, ii] = Quadcopter_tilted_arms_min_hover(ii,kf, m, g, Fdes, nmax, alpha0, n0, beta, theta, Display, Algorithm, maxIter);
                    nstar = round(dec*nstar)/dec;
                    alphastar = round(dec*alphastar)/dec;
                    % Test
                    if isequal(alphastar, alpha0) && isequal(nstar, n0)
                        worthH = worthH+1;
                    end
                    Hstar = m*g/(kf*norm(nstar)^2);
                    H0 = m*g/(kf*norm(n0)^2);
                    if Hstar < H0
                        Heff = [Heff H0]; % Hover efficiency in direction d
                        errH = errH+1;
                    else
                        Heff = [Heff Hstar]; % Hover efficiency in direction d
                    end
                else
                    alphastar = alpha0;
                    alphastar = round(dec*alphastar)/dec;
                    nstar = n0;
                    nstar = round(dec*nstar)/dec;
                    Hstar = m*g/(kf*norm(nstar)^2);
                    Heff = [Heff Hstar]; % Hover efficiency in direction d
                end
                % Test
%                 alphastar = round(dec*alphastar)/dec;
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta,nstar, L, g, Mb, Mp, R, false);
%                 pdotdot = round(dec*pdotdot)/dec;
%                 FH = [FH (m*pdotdot-m*g*d)];                
%                 staropt = [round(nstar); round(rad2deg(alphastar)).'; 0; round(100*(m*g/(kf*norm(nstar)^2)))];
%                 StarOpt = [StarOpt staropt];
                
            else
                Heff = [Heff 0];
            end
        end
    end
end
i0 = find(~vecnorm(D));
D(:,i0) = [];
F(:,i0) = [];
Feff(:,i0) = [];
M(:,i0) = [];
Meff(:,i0) = [];
Heff(:,i0) = [];

% % Test
% Fcheck(:,i0) = [];
% Star(:,i0) = [];
% FH(:,i0) = [];
% StarOpt(:,i0) = [];
% MM(:,i0) = [];

D_unit = D./vecnorm(D);
[D_unit,ia,ic] = unique(D_unit.', 'stable', 'rows');
D = D(:,ia.');
F = F(:,ia.');
Feff = Feff(:,ia.');
Meff = Meff(:,ia.');
Heff= Heff(:,ia.');
M = M(:,ia.');

% % Test
% Fcheck = Fcheck(:,ia.');
% FH = FH(:,ia.');
% Star = Star(:,ia.');
% StarOpt = StarOpt(:,ia.');
% MM = MM(:,ia.');


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

% % Test
% MM = [M;MM;Star];
end
