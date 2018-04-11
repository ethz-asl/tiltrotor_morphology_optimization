function [theta] = Quadcopter_tilted_arms_find_theta(beta)
%QUADCOPTER_TILTED_ARMS_FIND_THETA find arms theta angles
%   find theta to keep the angles between the arms equals for any beta
%% compute angle theta
fun = @(x) abs(x(1)) + abs(x(2)) + abs(x(3)) + abs(x(4)); %minimize the oposits angle
% Condition  Ax <= b  
A = []; 
b = [];          
% condition: Aeq.x = beq
Aeq = [];
beq = [];
% condition: lb <= x <= ub 
lb = [-pi/3 -pi/3 -pi/3 -pi/3].'; % lower bound
ub = [pi/3 pi/3 pi/3 pi/3].'; % upper bound
% initial guess:
x0 = [0 0 0 0];
options = optimoptions('fmincon', 'Display', 'off','Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',10000);
options=optimoptions(options,'MaxIter',10000);
theta = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, @(x) nonlincon(x, beta),  options);

end
function [c,ceq] = nonlincon(x, beta)
Op1 = rotz(rad2deg(x(1)))*roty(rad2deg(beta(1)))*[1 0 0].';
Op2 = rotz(rad2deg(pi/2+x(2)))*roty(rad2deg(beta(2)))*[1 0 0].';
Op3 = rotz(rad2deg(pi+x(3)))*roty(rad2deg(beta(3)))*[1 0 0].';
Op4 = rotz(rad2deg(3*pi/2+x(4)))*roty(rad2deg(beta(4)))*[1 0 0].';
c = [];
ceq(1) = atan2(norm(cross(Op1,Op2)), dot(Op1,Op2)) -atan2(norm(cross(Op2,Op3)), dot(Op2,Op3));
ceq(2) = atan2(norm(cross(Op2,Op3)), dot(Op2,Op3)) -atan2(norm(cross(Op3,Op4)), dot(Op3,Op4));
ceq(3) = atan2(norm(cross(Op3,Op4)), dot(Op3,Op4)) -atan2(norm(cross(Op4,Op1)), dot(Op4,Op1));

% ceq(1) = subspace(Op1,Op2) - subspace(Op2, Op3);
% ceq(2) = subspace(Op2,Op3) - subspace(Op3, Op4);
% ceq(3) = subspace(Op3,Op4) - subspace(Op4, Op1);

% ceq(1) = acos(Op1.'*Op2/(norm(Op1)*norm(Op2))) - acos(Op2.'*Op3/(norm(Op2)*norm(Op3)));
% ceq(2) = acos(Op2.'*Op3/(norm(Op2)*norm(Op3))) - acos(Op3.'*Op4/(norm(Op3)*norm(Op4)));
% ceq(3) = acos(Op3.'*Op4/(norm(Op3)*norm(Op4))) - acos(Op4.'*Op1/(norm(Op4)*norm(Op1)));
end

