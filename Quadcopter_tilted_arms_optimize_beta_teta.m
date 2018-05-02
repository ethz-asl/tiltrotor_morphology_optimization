function [betastar, tetastar, exitflag] = Quadcopter_tilted_arms_optimize_beta_teta(kf, km, L, g, Mb, Mp, R, nmin, nmax, betamin, betamax, thetamin, thetamax, alphadotmax, alphamin, alphamax, opt_iterations, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
%function [betastar, tetastar, exitflag] = Quadcopter_tilted_arms_optimize_beta_teta(kf, km, L, d, g, Mb, Mp, R, betamin, betamax, thetamin, thetamax, step, beta0, theta0, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance)
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
%% Fmincom optimisation using an other fmincom optimization.
beta = [x(1), x(2), x(3), x(3)];
theta = [x(4), x(5), x(6), x(7)];
[wRb, D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff, C, worthF, worthM, worthH, worthC, number_of_directions, TRI, F_surf, F_vol, M_surf, M_vol] = Quadcopter_tilted_arms_compute_metrics(beta ,theta, L, Mb, Mp, R, kf, km, nmin, nmax, alphamin, alphamax, g, step, false, Display, Algorithm, maxIter, StepTolerance, ConstraintTolerance,  opt_iterations, alphadotmax);

fun = -Fmin -Mmin - Hmin;

% %% Static matrix methode
% % The static matrix are static allocation matrix that do not depend on the
% % rotor orientation and speed.
% % (static matrix found using the file:Quadcopter_tilted_arms_Find_Static_Matrix.m)
% 
% % Vector containing all the decomposed vertical and horizontal forces:
% % Fdec = [kf*cos(alpha(1))*n(1)^2; kf*sin(alpha(1))*n(1)^2; kf*cos(alpha(2))*n(2)^2; 
% %         kf*sin(alpha(2))*n(2)^2; kf*cos(alpha(3))*n(3)^2; kf*sin(alpha(3))*n(3)^2; 
% %         kf*cos(alpha(4))*n(4)^2; kf*sin(alpha(4))*n(4)^2];
% 
% % This static matrix links Fdec to the force applied by the propellers to
% % the drone body 
% % F = m*p'' = A_F_static*Fdec
% A_F_static = [sin(x(1))*cos(x(5)), sin(x(5)), -sin(x(2))*sin(x(6)), ...
%               cos(x(6)), -sin(x(3))*cos(x(7)), -sin(x(7)), ...
%               sin(x(4))*sin(x(8)), -cos(x(8)); ...
%               sin(x(1))*sin(x(5)), -cos(x(5)), sin(x(2))*cos(x(6)), ...
%               sin(x(6)), -sin(x(3))*sin(x(7)), cos(x(7)), ...
%               -sin(x(4))*cos(x(8)), -sin(x(8)); ...
%               cos(x(1)), 0, cos(x(2)), 0,  cos(x(3)), 0, cos(x(4)), 0];
% 
% % This static matrix links Fdec to the torque applied by the propellers to
% % the drone body 
% % M = Ib*wb' = A_M_static*Fdec
% A_M_static =[L*sin(x(5))-km*sin(x(1))*cos(x(5))/kf, -L*sin(x(1))*cos(x(5))-km*sin(x(5))/kf, ...
%              L*cos(x(6))-km*sin(x(2))*sin(x(6))/kf, L*sin(x(2))*sin(x(6))+km*cos(x(6))/kf, ...
%              -L*sin(x(7))+km*sin(x(3))*cos(x(7))/kf, L*sin(x(3))*cos(x(7))+km*sin(x(7))/kf, ...
%              -L*cos(x(8))-km*sin(x(4))*sin(x(8))/kf, -L*sin(x(4))*sin(x(8))-km*cos(x(8))/kf; ...
%              -L*cos(x(5))-km*sin(x(1))*sin(x(5))/kf, -L*sin(x(1))*sin(x(5))+km*cos(x(5))/kf, ...
%              L*sin(x(6))+km*sin(x(2))*cos(x(6))/kf, -L*sin(x(2))*cos(x(6))+km*sin(x(6))/kf, ...
%              L*cos(x(7))+km*sin(x(3))*sin(x(7))/kf, L*sin(x(3))*sin(x(7))-km*cos(x(7))/kf, ...
%              -L*sin(x(8))-km*sin(x(4))*cos(x(8))/kf, L*sin(x(4))*cos(x(8))-km*sin(x(8))/kf; ...
%              -km*cos(x(1))/kf, -L*cos(x(1)), km*cos(x(2))/kf, -L*cos(x(2)), -km*cos(x(3))/kf, ...
%              -L*cos(x(3)), km*cos(x(4))/kf, -L*cos(x(4))];
% 
% % The Moore-Penrose pseudo inverse of the static matrices allow to find
% % Fdec from a desired force or torque applied on the drone.
% % The rotor orientation and speed can then be deduced from Fdec
% 
% % => Fdec = inv(A_F_static)*Fdes
% A_F_staticinv = pinv(A_F_static);
% 
% % => Fdec = inv(A_M_static)*Mdes
% A_M_staticinv = pinv(A_M_static);
% 
% %% Loop to create the direction matrix D with the desired number of direction:
% D = []; % Matrix containing every direction for which we want to compute the metrics
% F = []; % Matrix containing the maximum force appliable by the design in every direction of D
% M = []; % Vector containing maximum torque appliable by the design in every direction of D
% alphaF = [];
% nF = []; 
% alphaM = [];
% nM = []; 
% for i = 1:-step:-1
%     for j = 1:-step:-1
%         for k = 1:-step:-1% number of directions depends only on step size
%             d = [i j k].';
%             D = [D d];
%         end
%     end
% end
% i0 = find(~vecnorm(D)); 
% D(:,i0) = [];% Eliminate the [0; 0; 0] direction
% D_unit = D./vecnorm(D); % Create a normalized matrix of direction.
% D_unit = round(D_unit*10^5)/10^5;
% [D_unit,ia,ic] = unique(D_unit.', 'stable', 'rows');
% D_unit = D_unit.'; % Eliminate redundant directions in normalized D
% D = D(:,ia.');% Eliminate redundant directions in D
% for d = D_unit
%     %% find the max force in direction d using static matrix
%     Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers in direction d
%     
%     Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
%     
%     % Inverse substitution :
%     %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
%     %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
%     nF0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
%     nF0 = nF0/vecnorm(nF0);
%     nF0 = nmax^2*nF0/max(nF0); % nstar <= nmax
%     nF0(nF0<nmin^2) = nmin^2; % nmin <= nstar
%     nF0 = sqrt(nF0);
%     alphaF0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
%     % calculate angular and linear acceleration with this alphastar and nstar
%     [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphaF0, [x(1), x(2), x(3), x(4)], [x(5), x(6), x(7), x(8)],nF0, L, g, Mb, Mp, R, false);
%     F0 = m*pdotdot;
%     FN0 = norm(F0);
%     F = [F FN0];
%     alphaF = [alphaF; alphaF0];
%     nF = [nF; nF0.']; 
%     
%     %% find max torque in direction d
%     % find initial alpha and n for the optimisation to find the max torque in direction d
%     Mdes = d*(4*L*nmax^2*kf); % set desired torque to be equal to the maximal torque of the four propellers
%     Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
%     % Inverse substitution :
%     %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
%     %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
%     nM0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
%     nM0 = nM0/vecnorm(nM0);
%     nM0 = nmax^2*nM0/max(nM0); % 0 <= nstar <= nmax
%     nM0(nM0<nmin^2) = nmin^2;
%     nM0 = sqrt(nM0);
%     alphaM0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
%     % calculate angular and linear acceleration with this alphastar and nstar
%     [m, Ib, pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphaM0, [x(1), x(2), x(3), x(4)], [x(5), x(6), x(7), x(8)], nM0, L, g, Mb, Mp, R, false);
%     M0 = Ib*wbdot;
%     MN0 = norm(M0);
%     M = [M MN0];
%     alphaM = [alphaM; alphaM0];
%     nM = [nM; nM0.']; 
% end
% IFmin = find(F == min(F));
% IFmax = find(F == max(F));
% IMmin = find(M == min(M));
% IMmax = find(M == max(M));
% 
% alphaFmin = alphaF(IFmin(1),:);
% alphaFmax = alphaF(IFmax(1),:);
% alphaMmin = alphaM(IMmin(1),:);
% alphaMmax = alphaM(IMmax(1),:);
% 
% nFmin = nF(IFmin(1),:);
% nFmax = nF(IFmax(1),:);
% nMmin = nM(IMmin(1),:);
% nMmax = nM(IMmax(1),:);
% 
% paramFmin = [ alphaFmin, nFmin];
% paramFmax = [ alphaFmax, nFmax];
% paramMmin = [ alphaMmin, nMmin];
% paramMmax = [ alphaMmax, nMmax];
% 
% % Objective function (maximize the 2-norm of the torque produced by the propellers in the body frame):
% fun = -sqrt((-L*kf*paramMmin(5)^2*(sin(x(1))*sin(paramMmin(1))*cos(x(5)) - cos(paramMmin(1))*sin(x(5))) ...
%                   +L*kf*paramMmin(6)^2*(sin(x(2))*sin(paramMmin(2))* sin(x(6)) + cos(paramMmin(2))*cos(x(6))) ...
%                   +L*kf*paramMmin(7)^2*(sin(x(3))*sin(paramMmin(3))*cos(x(7)) - cos(paramMmin(3))* sin(x(7))) ...
%                   -L*kf*paramMmin(8)^2*(sin(x(4))*sin(paramMmin(4))*sin(x(8)) + cos(paramMmin(4))*cos(x(8))) ...
%                   -km*(sin(paramMmin(1))*sin(x(5)) + cos(paramMmin(1))*sin(x(1))*cos(x(5)))*paramMmin(5)^2 ...
%                   +km*(sin(paramMmin(2))*(cos(x(6))) - cos(paramMmin(2))*sin(x(2))*sin(x(6)))*paramMmin(6)^2 ...
%                   +km*(sin(paramMmin(3))*sin(x(7)) + cos(paramMmin(3))*sin(x(3))*cos(x(7)))*paramMmin(7)^2 ...
%                   -km*(sin(paramMmin(4))*(cos(x(8))) - cos(paramMmin(4))*sin(x(4))*sin(x(8)))*paramMmin(8)^2)^2 ...
%                  +(-L*kf*paramMmin(5)^2*(sin(x(1))*sin(paramMmin(1))*sin(x(5)) + cos(paramMmin(1))*cos(x(5))) ...
%                    -L*kf*paramMmin(6)^2*(sin(x(2))*sin(paramMmin(2))*cos(x(6)) - cos(paramMmin(2))* sin(x(6))) ...
%                    +L*kf*paramMmin(7)^2*(sin(x(3))*sin(paramMmin(3))* sin(x(7)) + cos(paramMmin(3))*cos(x(7))) ...
%                    +L*kf*paramMmin(8)^2*(sin(x(4))*sin(paramMmin(4))*cos(x(8)) - cos(paramMmin(4))*sin(x(8))) ...
%                    +km*(sin(paramMmin(1))*cos(x(5)) - cos(paramMmin(1))*sin(x(1))*sin(x(5)))*paramMmin(5)^2 ...
%                    +km*(sin(paramMmin(2))*sin(x(6)) + cos(paramMmin(2))*sin(x(2))*(cos(x(6))))*paramMmin(6)^2 ...
%                    -km*(sin(paramMmin(3))*cos(x(7)) - cos(paramMmin(3))*sin(x(3))*sin(x(7)))*paramMmin(7)^2 ...
%                    -km*(sin(paramMmin(4))*sin(x(8)) + cos(paramMmin(4))*sin(x(4))*cos(x(8)))*paramMmin(8)^2)^2 ...
%                  +(-L*kf*paramMmin(5)^2*cos(x(1))*sin(paramMmin(1)) - L*kf*paramMmin(6)^2*cos(x(2))*sin(paramMmin(2)) ...
%                    -L*kf*paramMmin(7)^2*cos(x(3))*sin(paramMmin(3)) - L*kf*paramMmin(8)^2*cos(x(4))*sin(paramMmin(4))...
%                    -km*cos(paramMmin(1))*cos(x(1))*paramMmin(5)^2 + km*cos(paramMmin(2))*cos(x(2))*paramMmin(6)^2 ...
%                    -km*cos(paramMmin(3))*cos(x(3))*paramMmin(7)^2 + km*cos(paramMmin(4))*cos(x(4))*paramMmin(8)^2)^2) ...
%            -sqrt((kf*(sin(paramFmin(1))*sin(x(5)) + cos(paramFmin(1))*sin(x(1))*cos(x(5)))*paramFmin(5)^2 ...
%                   +kf*(sin(paramFmin(2))*cos(x(6)) - cos(paramFmin(2))*sin(x(2))*sin(x(6)))*paramFmin(6)^2 ...
%                   -kf*(sin(paramFmin(3))*sin(x(7)) + cos(paramFmin(3))*sin(x(3))*cos(x(7)))*paramFmin(7)^2 ...
%                   -kf*(sin(paramFmin(4))*cos(x(8)) - cos(paramFmin(4))*sin(x(4))*sin(x(8)))*paramFmin(8)^2)^2 ...
%                  +(-kf*(sin(paramFmin(1))*cos(x(5)) - cos(paramFmin(1))*sin(x(1))*sin(x(5)))*paramFmin(5)^2 ...
%                    +kf*(sin(paramFmin(2))*sin(x(6)) + cos(paramFmin(2))*sin(x(2))*cos(x(6)))*paramFmin(6)^2 ...
%                    +kf*(sin(paramFmin(3))*(cos(x(7))) - cos(paramFmin(3))*sin(x(3))*sin(x(7)))*paramFmin(7)^2 ...
%                    -kf*(sin(paramFmin(4))*sin(x(8)) + cos(paramFmin(4))*sin(x(4))*cos(x(8)))*paramFmin(8)^2)^2 ...
%                  +(kf*cos(paramFmin(1))*cos(x(1))*paramFmin(5)^2 + kf*cos(paramFmin(2))*cos(x(2))*paramFmin(6)^2 ...
%                    +kf*cos(paramFmin(3))*cos(x(3))*paramFmin(7)^2+ kf*cos(paramFmin(4))*cos(x(4))*paramFmin(8)^2)^2)...
%            -sqrt((-L*kf*paramMmax(5)^2*(sin(x(1))*sin(paramMmax(1))*cos(x(5)) - cos(paramMmax(1))*sin(x(5))) ...
%                   +L*kf*paramMmax(6)^2*(sin(x(2))*sin(paramMmax(2))* sin(x(6)) + cos(paramMmax(2))*cos(x(6))) ...
%                   +L*kf*paramMmax(7)^2*(sin(x(3))*sin(paramMmax(3))*cos(x(7)) - cos(paramMmax(3))* sin(x(7))) ...
%                   -L*kf*paramMmax(8)^2*(sin(x(4))*sin(paramMmax(4))*sin(x(8)) + cos(paramMmax(4))*cos(x(8))) ...
%                   -km*(sin(paramMmax(1))*sin(x(5)) + cos(paramMmax(1))*sin(x(1))*cos(x(5)))*paramMmax(5)^2 ...
%                   +km*(sin(paramMmax(2))*(cos(x(6))) - cos(paramMmax(2))*sin(x(2))*sin(x(6)))*paramMmax(6)^2 ...
%                   +km*(sin(paramMmax(3))*sin(x(7)) + cos(paramMmax(3))*sin(x(3))*cos(x(7)))*paramMmax(7)^2 ...
%                   -km*(sin(paramMmax(4))*(cos(x(8))) - cos(paramMmax(4))*sin(x(4))*sin(x(8)))*paramMmax(8)^2)^2 ...
%                  +(-L*kf*paramMmax(5)^2*(sin(x(1))*sin(paramMmax(1))*sin(x(5)) + cos(paramMmax(1))*cos(x(5))) ...
%                    -L*kf*paramMmax(6)^2*(sin(x(2))*sin(paramMmax(2))*cos(x(6)) - cos(paramMmax(2))* sin(x(6))) ...
%                    +L*kf*paramMmax(7)^2*(sin(x(3))*sin(paramMmax(3))* sin(x(7)) + cos(paramMmax(3))*cos(x(7))) ...
%                    +L*kf*paramMmax(8)^2*(sin(x(4))*sin(paramMmax(4))*cos(x(8)) - cos(paramMmax(4))*sin(x(8))) ...
%                    +km*(sin(paramMmax(1))*cos(x(5)) - cos(paramMmax(1))*sin(x(1))*sin(x(5)))*paramMmax(5)^2 ...
%                    +km*(sin(paramMmax(2))*sin(x(6)) + cos(paramMmax(2))*sin(x(2))*(cos(x(6))))*paramMmax(6)^2 ...
%                    -km*(sin(paramMmax(3))*cos(x(7)) - cos(paramMmax(3))*sin(x(3))*sin(x(7)))*paramMmax(7)^2 ...
%                    -km*(sin(paramMmax(4))*sin(x(8)) + cos(paramMmax(4))*sin(x(4))*cos(x(8)))*paramMmax(8)^2)^2 ...
%                  +(-L*kf*paramMmax(5)^2*cos(x(1))*sin(paramMmax(1)) - L*kf*paramMmax(6)^2*cos(x(2))*sin(paramMmax(2)) ...
%                    -L*kf*paramMmax(7)^2*cos(x(3))*sin(paramMmax(3)) - L*kf*paramMmax(8)^2*cos(x(4))*sin(paramMmax(4))...
%                    -km*cos(paramMmax(1))*cos(x(1))*paramMmax(5)^2 + km*cos(paramMmax(2))*cos(x(2))*paramMmax(6)^2 ...
%                    -km*cos(paramMmax(3))*cos(x(3))*paramMmax(7)^2 + km*cos(paramMmax(4))*cos(x(4))*paramMmax(8)^2)^2) ...
%            -sqrt((kf*(sin(paramFmax(1))*sin(x(5)) + cos(paramFmax(1))*sin(x(1))*cos(x(5)))*paramFmax(5)^2 ...
%                   +kf*(sin(paramFmax(2))*cos(x(6)) - cos(paramFmax(2))*sin(x(2))*sin(x(6)))*paramFmax(6)^2 ...
%                   -kf*(sin(paramFmax(3))*sin(x(7)) + cos(paramFmax(3))*sin(x(3))*cos(x(7)))*paramFmax(7)^2 ...
%                   -kf*(sin(paramFmax(4))*cos(x(8)) - cos(paramFmax(4))*sin(x(4))*sin(x(8)))*paramFmax(8)^2)^2 ...
%                  +(-kf*(sin(paramFmax(1))*cos(x(5)) - cos(paramFmax(1))*sin(x(1))*sin(x(5)))*paramFmax(5)^2 ...
%                    +kf*(sin(paramFmax(2))*sin(x(6)) + cos(paramFmax(2))*sin(x(2))*cos(x(6)))*paramFmax(6)^2 ...
%                    +kf*(sin(paramFmax(3))*(cos(x(7))) - cos(paramFmax(3))*sin(x(3))*sin(x(7)))*paramFmax(7)^2 ...
%                    -kf*(sin(paramFmax(4))*sin(x(8)) + cos(paramFmax(4))*sin(x(4))*cos(x(8)))*paramFmax(8)^2)^2 ...
%                  +(kf*cos(paramFmax(1))*cos(x(1))*paramFmax(5)^2 + kf*cos(paramFmax(2))*cos(x(2))*paramFmax(6)^2 ...
%                    +kf*cos(paramFmax(3))*cos(x(3))*paramFmax(7)^2+ kf*cos(paramFmax(4))*cos(x(4))*paramFmax(8)^2)^2);
end
end