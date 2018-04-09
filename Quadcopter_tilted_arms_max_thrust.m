function [alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, nhover, alpha0, n0, d, beta, theta)
%[alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, nhover, alpha0, n0, d, beta, theta)
%QUADCOPTER_TILTED_ARMS_MAX_THRUST find optimal tilting angles and rotor
%speed

%% Optimization of alpha and n
% maximize norm of the Thrust F in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]                                                                                                                                                               
% Condition  Ax <= b        
A = []; 
b = [];                 

% condition: Aeq.x = beq     
Aeq = [];
beq = [];

% condition: lb <= x <= ub 
lb = [-pi -pi -pi -pi 0 0 0 0].'; % lower bound
ub = [pi pi pi pi nmax nmax nmax nmax].'; % upper bound

% initial guess:
x0 = [alpha0, n0.'];

options = optimoptions('fmincon', 'Display', 'off','Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
x = fmincon(@(x) ThrustNorm(kf, theta, beta, x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf, theta, beta, d),  options);

% Substitution:
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
k = 10.^Ndecimals;
nstar = round(k*nstar)/k;
alphastar = round(k*alphastar)/k;
end
function [squareNormT] = ThrustNorm(kf, beta, theta, x)
T = Thrust(kf, theta, beta, x);
% Objecive function maximize the squared norm of the thrust T:
squareNormT = -(T(1)^2+ T(2)^2+ T(3)^2);
end
function [T] = Thrust(kf, theta, beta, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% => T = [kf*(sin(alpha1)*sin(theta1) + cos(alpha1)*sin(beta1)*cos(theta1))*n1^2 + kf*(sin(alpha2)*sin(theta2 + pi/2) + cos(alpha2)*sin(beta2)*cos(theta2 + pi/2))*n2^2 - kf*(sin(alpha3)*sin(theta3) + cos(alpha3)*sin(beta3)*cos(theta3))*n3^2 + kf*(sin(alpha4)*sin(theta4 + (3*pi)/2) + cos(alpha4)*sin(beta4)*cos(theta4 + (3*pi)/2))*n4^2
%         -kf*(sin(alpha1)*cos(theta1) - cos(alpha1)*sin(beta1)*sin(theta1))*n1^2 - kf*(sin(alpha2)*cos(theta2 + pi/2) - cos(alpha2)*sin(beta2)*sin(theta2 + pi/2))*n2^2 + kf*(sin(alpha3)*cos(theta3) - cos(alpha3)*sin(beta3)*sin(theta3))*n3^2 - kf*(sin(alpha4)*cos(theta4 + (3*pi)/2) - cos(alpha4)*sin(beta4)*sin(theta4 + (3*pi)/2))*n4^2
%         kf*cos(alpha1)*cos(beta1)*n1^2 + kf*cos(alpha2)*cos(beta2)*n2^2 + kf*cos(alpha3)*cos(beta3)*n3^2 + kf*cos(alpha4)*cos(beta4)*n4^2] (Thrust T)
T1 = kf*(sin(x(1))*sin(theta(1)) + cos(x(1))*sin(beta(1))*cos(theta(1)))*x(5)^2 + kf*(sin(x(2))*sin(theta(2) + pi/2) + cos(x(2))*sin(beta(2))*cos(theta(2) + pi/2))*x(6)^2 - kf*(sin(x(3))*sin(theta(3)) + cos(x(3))*sin(beta(3))*cos(theta(3)))*x(7)^2 + kf*(sin(x(4))*sin(theta(4) + (3*pi)/2) + cos(x(4))*sin(beta(4))*cos(theta(4) + (3*pi)/2))*x(8)^2;
T2 = -kf*(sin(x(1))*cos(theta(1)) - cos(x(1))*sin(beta(1))*sin(theta(1)))*x(5)^2 - kf*(sin(x(2))*cos(theta(2) + pi/2) - cos(x(2))*sin(beta(2))*sin(theta(2) + pi/2))*x(6)^2 + kf*(sin(x(3))*cos(theta(3)) - cos(x(3))*sin(beta(3))*sin(theta(3)))*x(7)^2 - kf*(sin(x(4))*cos(theta(4) + (3*pi)/2) - cos(x(4))*sin(beta(4))*sin(theta(4) + (3*pi)/2))*x(8)^2;
T3 = kf*cos(x(1))*cos(beta(1))*x(5)^2 + kf*cos(x(2))*cos(beta(2))*x(6)^2 + kf*cos(x(3))*cos(beta(3))*x(7)^2 + kf*cos(x(4))*cos(beta(4))*x(8)^2;
T = [T1 T2 T3];
end
function [c,ceq] = nonlinconF(x, kf, theta, beta, d)
% function [c,ceq] = nonlincon(x, kf,d)
% % Hover condition
% if d(3) < 0
%     c = [];
% else
%     c(1) = -cos(x(1))*x(5)^2 - cos(x(2))*x(6)^2 - cos(x(3))*x(7)^2 - cos(x(4))*x(8)^2;
% end
% No hover cdt
c = [];
% Thrust parallel to d conditions: Ceq(x) = 0
T = Thrust(kf, theta, beta, x);
ceq(1) = T(2)*d(3) - T(3)*d(2);
ceq(2) = T(3)*d(1) - T(1)*d(3);
ceq(3) = T(1)*d(2) - T(2)*d(1);
end
