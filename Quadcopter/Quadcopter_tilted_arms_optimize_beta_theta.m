function [betastar, tetastar, exitflag] = Quadcopter_tilted_arms_optimize_beta_theta(kf, km, L, g, Mb, Mp, R, nmin, nmax, betamin, betamax, thetamin, thetamax, alphadotmax, alphamin, alphamax, opt_iterations, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%function [betastar, tetastar, exitflag] = Quadcopter_tilted_arms_optimize_beta_theta(kf, km, L, d, g, Mb, Mp, R, betamin, betamax, thetamin, thetamax, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%QUADCOPTER_TILTED_ARMS_MAX_TORQUE find optimal tilting angles and Rotor 
%
%speed  
%   Optimize alpha and n so the drone produce the maximal torque in arbitrary direction d

%% Optimization of alpha and n
% maximize norm of the Torque M in an arbitrairy direction d:

% x = [beta1 beta2 beta3 beta4 teta1 teta2 teta3 teta4]

% Condition  Ax <= b        
A = []; 
b = [];                 

% Condition: Aeq.x = beq     
Aeq = [];
beq = [];

% Condition: lb <= x <= ub 
lb = [betamin betamin betamin betamin thetamin thetamin thetamin thetamin].'; % lower bound
ub = [betamax betamax betamax betamax thetamax thetamax thetamax thetamax].'; % upper bound

% initial guess:
x0 = [beta0, theta0];
         
% optimization options       
options = optimoptions('fmincon', 'Display', Display, 'Algorithm',Algorithm, 'StepTolerance', StepTolerance, 'ConstraintTolerance', ConstraintTolerance);
options=optimoptions(options, 'MaxFunEvals', maxIter);
options=optimoptions(options,'MaxIter', maxIter);

% actual optimization
[xstar,fval,exitflag,output] = fmincon(@ (x) objective_function(x, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance,  opt_iterations, alphadotmax), x0, A, b, Aeq, beq, lb, ub,[],  options);

% Solution of the optimization
betastar = [xstar(1) xstar(2) xstar(3) xstar(4)];
tetastar = [xstar(5) xstar(6) xstar(7) xstar(8)];

%% Objective function function  
function [fun] = objective_function(x, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance,  opt_iterations, alphadotmax)
%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb0 = rotz(rad2deg(roll0))*roty(rad2deg(pitch0))*rotz(rad2deg(yaw0)); % Rotation Matrix mapping body frame to inertial frame
%% Fmincom optimisation using the static matrix

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
A_F_static = [sin(x(1))*cos(x(5)), sin(x(5)), -sin(x(2))*sin(x(6)), ...
              cos(x(6)), -sin(x(3))*cos(x(7)), -sin(x(7)), ...
              sin(x(4))*sin(x(8)), -cos(x(8)); ...
              sin(x(1))*sin(x(5)), -cos(x(5)), sin(x(2))*cos(x(6)), ...
              sin(x(6)), -sin(x(3))*sin(x(7)), cos(x(7)), ...
              -sin(x(4))*cos(x(8)), -sin(x(8)); ...
              cos(x(1)), 0, cos(x(2)), 0,  cos(x(3)), 0, cos(x(4)), 0];

% This static matrix links Fdec to the torque applied by the propellers to
% the drone body 
% M = Ib*wb' = A_M_static*Fdec
A_M_static =[L*sin(x(5))-km*sin(x(1))*cos(x(5))/kf, -L*sin(x(1))*cos(x(5))-km*sin(x(5))/kf, ...
             L*cos(x(6))-km*sin(x(2))*sin(x(6))/kf, L*sin(x(2))*sin(x(6))+km*cos(x(6))/kf, ...
             -L*sin(x(7))+km*sin(x(3))*cos(x(7))/kf, L*sin(x(3))*cos(x(7))+km*sin(x(7))/kf, ...
             -L*cos(x(8))-km*sin(x(4))*sin(x(8))/kf, -L*sin(x(4))*sin(x(8))-km*cos(x(8))/kf; ...
             -L*cos(x(5))-km*sin(x(1))*sin(x(5))/kf, -L*sin(x(1))*sin(x(5))+km*cos(x(5))/kf, ...
             L*sin(x(6))+km*sin(x(2))*cos(x(6))/kf, -L*sin(x(2))*cos(x(6))+km*sin(x(6))/kf, ...
             L*cos(x(7))+km*sin(x(3))*sin(x(7))/kf, L*sin(x(3))*sin(x(7))-km*cos(x(7))/kf, ...
             -L*sin(x(8))-km*sin(x(4))*cos(x(8))/kf, L*sin(x(4))*cos(x(8))-km*sin(x(8))/kf; ...
             -km*cos(x(1))/kf, -L*cos(x(1)), km*cos(x(2))/kf, -L*cos(x(2)), -km*cos(x(3))/kf, ...
             -L*cos(x(3)), km*cos(x(4))/kf, -L*cos(x(4))];

% The Moore-Penrose pseudo inverse of the static matrices allow to find
% Fdec from a desired force or torque applied on the drone.
% The rotor orientation and speed can then be deduced from Fdec

% => Fdec = inv(A_F_static)*Fdes
A_F_staticinv = pinv(A_F_static);

% => Fdec = inv(A_M_static)*Mdes
A_M_staticinv = pinv(A_M_static);

%% Loop to create the direction matrix D with the desired number of direction:
D = []; % Matrix containing every direction for which we want to compute the metrics
F = []; % Matrix containing the maximum force appliable by the design in every direction of D
M = []; % Vector containing maximum torque appliable by the design in every direction of D
Heff = [];
alphaF = [];
nF = []; 
alphaM = [];
nM = []; 
alphaH = [];
nH = []; 
for i = 1:-step:-1
    for j = 1:-step:-1
        for k = 1:-step:-1% number of directions depends only on step size
            d = [i j k].';
            D = [D d];
        end
    end
end
i0 = find(~vecnorm(D)); 
D(:,i0) = [];% Eliminate the [0; 0; 0] direction
D_unit = D./vecnorm(D); % Create a normalized matrix of direction.
D_unit = round(D_unit*10^5)/10^5;
[D_unit,ia,ic] = unique(D_unit.', 'stable', 'rows');
D_unit = D_unit.'; % Eliminate redundant directions in normalized D
D = D(:,ia.');% Eliminate redundant directions in D
for d = D_unit
    %% find the max force in direction d using static matrix
    Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers in direction d
    
    Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
    
    % Inverse substitution :
    %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    nF0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
    nF0 = nF0/vecnorm(nF0);
    nF0 = nmax^2*nF0/max(nF0); % nstar <= nmax
    nF0(nF0<nmin^2) = nmin^2; % nmin <= nstar
    nF0 = sqrt(nF0);
    alphaF0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
    % calculate angular and linear acceleration with this alphastar and nstar
    [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphaF0, [x(1), x(2), x(3), x(4)], [x(5), x(6), x(7), x(8)],nF0, L, g, Mb, Mp, R, false);
    F0 = m*pdotdot;
    FN0 = norm(F0);
    F = [F FN0];
    alphaF = [alphaF; alphaF0];
    nF = [nF; nF0.']; 
    
    %% find max torque in direction d
    % find initial alpha and n for the optimisation to find the max torque in direction d
    Mdes = d*(4*L*nmax^2*kf); % set desired torque to be equal to the maximal torque of the four propellers
    Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
    % Inverse substitution :
    %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    nM0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
    nM0 = nM0/vecnorm(nM0);
    nM0 = nmax^2*nM0/max(nM0); % 0 <= nstar <= nmax
    nM0(nM0<nmin^2) = nmin^2;
    nM0 = sqrt(nM0);
    alphaM0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
    % calculate angular and linear acceleration with this alphastar and nstar
    [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphaM0, [x(1), x(2), x(3), x(4)], [x(5), x(6), x(7), x(8)], nM0, L, g, Mb, Mp, R, false);
    M0 = Ib*wbdot;
    MN0 = norm(M0);
    M = [M MN0];
    alphaM = [alphaM; alphaM0];
    nM = [nM; nM0.'];
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
    H0 = m*g/(kf*norm(n0)^2);
    Heff = [Heff H0]; % Hover efficiency in direction d
    alphaH = [alphaH; alpha0];
    nH = [nH; n0.'];
end
Heff = 100*Heff;
IFmin = find(F == min(F));
IFmax = find(F == max(F));
IMmin = find(M == min(M));
IMmax = find(M == max(M));
IHmin = find(Heff == min(Heff));
IHmax = find(Heff == max(Heff));

alphaFmin = alphaF(IFmin(1),:);
alphaFmax = alphaF(IFmax(1),:);
alphaMmin = alphaM(IMmin(1),:);
alphaMmax = alphaM(IMmax(1),:);

nFmin = nF(IFmin(1),:);
nFmax = nF(IFmax(1),:);
nMmin = nM(IMmin(1),:);
nMmax = nM(IMmax(1),:);
nHmin = nM(IHmin(1),:);
nHmax = nM(IHmax(1),:);

% Objective function (maximize the 2-norm of the minimum guaranteed torque produced by the propellers in the body frame):
fun = -4*norm(M(IMmin(1))) -20*norm(F(IFmin(1)))  -norm(Heff(IHmin(1)));

end
% function [M] = get_torque(L, kf, km, beta, theta, alpha, n)
%     % Torque produced by the propellers in the body frame
%     M = [(- L*kf*n(1)^2*(sin(beta(1))*sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(theta(1))) ...
%         + L*kf*n(2)^2*(sin(beta(2))*sin(alpha(2))* sin(theta(2)) + cos(alpha(2))*cos(theta(2))) ...
%         + L*kf*n(3)^2*(sin(beta(3))*sin(alpha(3))*cos(theta(3)) - cos(alpha(3))* sin(theta(3))) ...
%         - L*kf*n(4)^2*(sin(beta(4))*sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*cos(theta(4))) ...
%         -km*(sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*sin(beta(1))*cos(theta(1)))*n(1)^2 ...
%         + km*(sin(alpha(2))*(cos(theta(2))) - cos(alpha(2))*sin(beta(2))*sin(theta(2)))*n(2)^2 ...
%         + km*(sin(alpha(3))*sin(theta(3)) + cos(alpha(3))*sin(beta(3))*cos(theta(3)))*n(3)^2 ...
%         - km*(sin(alpha(4))*(cos(theta(4))) - cos(alpha(4))*sin(beta(4))*sin(theta(4)))*n(4)^2); ...
%         (- L*kf*n(1)^2*(sin(beta(1))*sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*cos(theta(1))) ...
%         - L*kf*n(2)^2*(sin(beta(2))*sin(alpha(2))*cos(theta(2)) - cos(alpha(2))* sin(theta(2))) ...
%         + L*kf*n(3)^2*(sin(beta(3))*sin(alpha(3))* sin(theta(3)) + cos(alpha(3))*cos(theta(3))) ...
%         + L*kf*n(4)^2*(sin(beta(4))*sin(alpha(4))*cos(theta(4)) - cos(alpha(4))*sin(theta(4))) ...
%         + km*(sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(beta(1))*sin(theta(1)))*n(1)^2 ...
%         + km*(sin(alpha(2))*sin(theta(2)) + cos(alpha(2))*sin(beta(2))*(cos(theta(2))))*n(2)^2 ...
%         - km*(sin(alpha(3))*cos(theta(3)) - cos(alpha(3))*sin(beta(3))*sin(theta(3)))*n(3)^2 ...
%         - km*(sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*sin(beta(4))*cos(theta(4)))*n(4)^2); ...
%         (- L*kf*n(1)^2*cos(beta(1))*sin(alpha(1)) - L*kf*n(2)^2*cos(beta(2))*sin(alpha(2)) ...
%         - L*kf*n(3)^2*cos(beta(3))*sin(alpha(3)) - L*kf*n(4)^2*cos(beta(4))*sin(alpha(4))...
%         -km*cos(alpha(1))*cos(beta(1))*n(1)^2 + km*cos(alpha(2))*cos(beta(2))*n(2)^2 ...
%         - km*cos(alpha(3))*cos(beta(3))*n(3)^2 + km*cos(alpha(4))*cos(beta(4))*n(4)^2)];
% end
% function [T] = get_thrust(kf, beta, theta, alpha, n)
%         % Thrusts produced by the propellers in the body frame
%         T = [(kf*(sin(alpha(1))*sin(theta(1)) + cos(alpha(1))*sin(beta(1))*cos(theta(1)))*n(1)^2 ...
%             +kf*(sin(alpha(2))*cos(theta(2)) - cos(alpha(2))*sin(beta(2))*sin(theta(2)))*n(2)^2 ...
%             -kf*(sin(alpha(3))*sin(theta(3)) + cos(alpha(3))*sin(beta(3))*cos(theta(3)))*n(3)^2 ...
%             -kf*(sin(alpha(4))*cos(theta(4)) - cos(alpha(4))*sin(beta(4))*sin(theta(4)))*n(4)^2); ...
%             (-kf*(sin(alpha(1))*cos(theta(1)) - cos(alpha(1))*sin(beta(1))*sin(theta(1)))*n(1)^2 ...
%             +kf*(sin(alpha(2))*sin(theta(2)) + cos(alpha(2))*sin(beta(2))*cos(theta(2)))*n(2)^2 ...
%             +kf*(sin(alpha(3))*(cos(theta(3))) - cos(alpha(3))*sin(beta(3))*sin(theta(3)))*n(3)^2 ...
%             -kf*(sin(alpha(4))*sin(theta(4)) + cos(alpha(4))*sin(beta(4))*cos(theta(4)))*n(4)^2); ...
%             (kf*cos(alpha(1))*cos(beta(1))*n(1)^2 + kf*cos(alpha(2))*cos(beta(2))*n(2)^2 ...
%             +kf*cos(alpha(3))*cos(beta(3))*n(3)^2+ kf*cos(alpha(4))*cos(beta(4))*n(4)^2)];
% end
end