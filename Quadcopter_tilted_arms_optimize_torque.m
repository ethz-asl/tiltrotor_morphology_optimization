function [alphastar, nstar] = Quadcopter_optimize_torque(kf, nmax, nhover, d, alpha0, n0)
%QUADCOPTER_OPTIMISE_TORQUE Summary of this function goes here
%   Detailed explanation goes here

%% Optimization of alpha and n
% maximize norm of the Torque M in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% => M = [(L*kf*cos(x(2)) - km*sin(x(2)))*x(6)^2 - (L*kf*cos(x(4)) - sin(x(4)))*x(8)^2; ...
%         -(L*kf*cos(x(1)) + km*sin(x(1)))*x(5)^2 + (L*kf*cos(x(3)) + sin(x(3)))*x(7)^2; ...
%         -(L*kf*sin(x(1)) - km*cos(x(1)))*x(5)^2 - (L*kf*sin(x(2)) + km*cos(x(2)))*x(6)^2 ...
%         -(L*kf*sin(x(3)) - km*cos(x(3)))*x(7)^2 - (L*kf*sin(x(4)) + km*cos(x(4)))*x(8)^2] (Thrust T)
% Neglect propeller drag:
fun = @(x)-sqrt((L*kf*cos(x(2))*x(6)^2 - L*kf*cos(x(4))*x(8)^2)^2 + ...
                 (-L*kf*cos(x(1))*x(5)^2 + L*kf*cos(x(3))*x(7)^2)^2 +...
                 (-L*kf*sin(x(1))*x(5)^2 - L*kf*sin(x(2))*x(6)^2 ...
                 -L*kf*sin(x(3))*x(7)^2 - L*kf*sin(x(4))*x(8)^2)^2);
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

options = optimoptions('fmincon','Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
% options = optimoptions(options, 'fmincon','Display','iter','Algorithm','sqp');
x = fmincon(@(x)fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconM(x, kf, km,d, L),  options);

% Substitution:
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
k = 10.^Ndecimals;
nstar = round(k*nstar)/k;
alphastar = round(k*alphastar)/k;


end

function [c,ceq] = nonlinconM(x, kf,km, d, L)
% function [c,ceq] = nonlincon(x, kf, nmax,d, L)

c = [];

% Neglect propeller drag:
%Torque
M1 = L*kf*cos(x(2))*x(6)^2 - L*kf*cos(x(4))*x(8)^2;
M2 = -L*kf*cos(x(1))*x(5)^2 + L*kf*cos(x(3))*x(7)^2; 
M3 = -L*kf*sin(x(1))*x(5)^2 - L*kf*sin(x(2))*x(6)^2 -L*kf*sin(x(3))*x(7)^2 - L*kf*sin(x(4))*x(8)^2;

% consider propeller drag:
% %Torque
% M1 = (L*kf*cos(x(2)) - km*sin(x(2)))*x(6)^2 - (L*kf*cos(x(4)) - sin(x(4)))*x(8)^2;
% M2 = -(L*kf*cos(x(1)) + km*sin(x(1)))*x(5)^2 + (L*kf*cos(x(3)) + sin(x(3)))*x(7)^2; 
% M3 = -(L*kf*sin(x(1)) - km*cos(x(1)))*x(5)^2 - (L*kf*sin(x(2)) + km*cos(x(2)))*x(6)^2 ...
%       -(L*kf*sin(x(3)) - km*cos(x(3)))*x(7)^2 - (L*kf*sin(x(4)) + km*cos(x(4)))*x(8)^2;

% Torque vector parallel to d condition:
ceq(1) = M2*d(3) - M3*d(2);
ceq(2) = M3*d(1) - M1*d(3);
ceq(3) = M1*d(2) - M2*d(1);
% Force applied on MAV equal to zero condition:
ceq(4) = sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2;
ceq(5) = -sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2;
ceq(6) = cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2;
end

