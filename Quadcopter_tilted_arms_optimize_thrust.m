function [alphastar, nstar] = Quadcopter_optimize_thrust(kf, nmax, nhover, alpha0, n0)
%QUADCOPTER_OPTIMIZE_THRUST Summary of this function goes here

%% Optimization of alpha and n
% maximize norm of the Thrust F in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% => T = kf*[sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2 ; -sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2; ...
%            cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2] (Thrust T)
fun = @(x)-sqrt((sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)^2 + (-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2)^2 ...
                 + (cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2);

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
x = fmincon(@(x)fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf,d),  options);

% Substitution:
alphastar = [x(1) x(2) x(3) x(4)];
nstar = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
k = 10.^Ndecimals;
nstar = round(k*nstar)/k;
alphastar = round(k*alphastar)/k;
end

function [c,ceq] = nonlinconF(x, kf,d)
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
ceq(1) = d(3)*(-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2)-d(2)*(cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2);
ceq(2) = -d(3)*(sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)+d(1)*(cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2);
ceq(3) = d(2)*(sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)-d(1)*(-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2);
end
