function [alphastar, nstar] = Quadcopter_tilted_arms_opt_hover(kf, f, nmax, alpha0, n0, beta, theta)
% [alphastar, nstar] = Quadcopter_tilted_arms_opt_hover(kf,m, g, nmax, alpha0, n0, beta, theta)
%QUADCOPTER_TILTED_ARMS_OPT_HOVER finds the optimal hover mode 
%   Optimize alpha and n to obtain the most efficient hover when drone is
%   oriented in direction d
%% Optimization of alpha and n
% maximize norm of the Torque M in an arbitrairy direction d:
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

fun = @(x) abs(x(5) + x(6) + x(7) + x(8));

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
[x,fval,exitflag,output] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf, theta, beta, f), options);
% Substitution:
if exitflag == 1
    alphastar = [x(1) x(2) x(3) x(4)];
    nstar = [x(5) x(6) x(7) x(8)].';
else
    alphastar = alpha0;
    nstar = n0;
end



end
function [T] = Thrust(kf, theta, beta, x)
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
alpha = [x(1); x(2); x(3); x(4)];
n = [x(5); x(6); x(7); x(8)];
bRp1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
bRp2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
bRp3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
bRp4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
Tp1 = [0 0 kf*n(1)^2].'; % Thrust vector propeller 1
Tp2 = [0 0 kf*n(2)^2].'; % Thrust vector propeller 2
Tp3 = [0 0 kf*n(3)^2].'; % Thrust vector propeller 3
Tp4 = [0 0 kf*n(4)^2].'; % Thrust vector propeller 4
T = (bRp1*Tp1 + bRp2*Tp2 + bRp3*Tp3 + bRp4*Tp4);
end
function [c,ceq] = nonlinconF(x, kf, theta, beta, f)
% function [c,ceq] = nonlincon(x, kf,d)
c = [];
% Thrust parallel to f conditions: Ceq(x) = 0
T = Thrust(kf, theta, beta, x);
ceq(1) = T(1)-f(1);
ceq(2) = T(2)-f(2);
ceq(3) = T(3)-f(3);
end


