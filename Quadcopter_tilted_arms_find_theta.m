function [theta] = Quadcopter_tilted_arms_find_theta(beta, optim)
%QUADCOPTER_TILTED_ARMS_FIND_THETA find arms theta angles
%   find theta to keep the angles between the arms equals for any beta
%% compute angle theta
if optim
    fun = @(x) abs(x(2) + x(3)); %minimize the oposits angle
    % Condition  Ax <= b        
    A = []; 
    b = [];                 
    % condition: Aeq.x = beq     
    Aeq = [];
    beq = [];
    % condition: lb <= x <= ub 
    lb = [-pi/2 -pi/2 -pi/2 -pi/2].'; % lower bound
    ub = [pi/2 pi/2 pi/2 pi/2].'; % upper bound
    % initial guess:
    x0 = [0 0 0 0];
    options = optimoptions('fmincon', 'Display', 'off','Algorithm','sqp');
    options=optimoptions(options, 'MaxFunEvals',100000);
    options=optimoptions(options,'MaxIter',100000);
    theta = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlincon(x, beta),  options);
else
    theta = [0 0 0 0];
end
end
function [c,ceq] = nonlincon(x, beta)
a1 = cos(beta(1))*cos(beta(2));
b1 = sin(beta(1))*sin(beta(2));
a2 = cos(beta(2))*cos(beta(3));
b2 = sin(beta(2))*sin(beta(3));
a3 = cos(beta(3))*cos(beta(4));
b3 = sin(beta(3))*sin(beta(4));
a4 = cos(beta(4))*cos(beta(1));
b4 = sin(beta(4))*sin(beta(1));
c = [];
% angles between arms are equals condition:
ceq(1) = acos(a1*sin(x(1)-x(2)) + b1) - acos(a2*sin(x(2)-x(3)) + b2);
ceq(2) = acos(a2*sin(x(2)-x(3)) + b2) - acos(a3*sin(x(3)-x(4)) + b3);
ceq(3) = acos(a3*sin(x(3)-x(4)) + b3) - acos(a4*sin(x(4)-x(1)) + b4);
end

