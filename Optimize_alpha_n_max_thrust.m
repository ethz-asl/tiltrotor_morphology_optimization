function [alpha, n] = Optimize_alpha_n_max_thrust(kf, nmax, nhover, g, d, m)
% [alpha, n] = Optimize_F_alpha_n(kf, nmax, nhover, g, d, m)
% OPTIMIZE_ALPHA_N_MAX_THRUST Optimization of alpha and n
%   maximize norm of the Thrust F in an arbitrairy direction d:

%% Optimization of alpha and n
% maximize norm of the Thrust F in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% => F = kf*[sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2 ; -sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2; ...
%            cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2 -m*g/kf] (Thrust T)
fun = @(x)-sqrt((sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)^2 + (-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2)^2 ...
                 + (cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2 -m*g/kf)^2);

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
if d(1) == 0 
    if d(2) == 0
        if d(3) == 0 
            % case: d = [0 0 0]
            fprintf(2,"Error, no direction selected")
            alpha = [0 0 0 0];
            n = [0 0 0 0].';
            return
        elseif d(3) < 0
            % case: d = [0 0 -z]
            % initial guess: max thrust -z direction:
            x0 = [pi pi pi pi nmax nmax nmax nmax];
        else 
            % case: d = [0 0 z]
            % initial guess: max thrust z direction:
            x0 = [0 0 0 0 nmax nmax nmax nmax];
        end
    elseif d(2) < 0
        if d(3) == 0 
            % case: d = [0 -y 0]
            % initial guess: max thrust y- direction:
            x0 = [pi/2 0 -pi/2 0 nmax 133.3593 nmax 133.3593];
        elseif d(3) < 0
            % case: d = [0 -y -z]
            % init guess is thrust -y -z direction:
            x0 = [pi/2 pi -pi/2 pi nmax nmax nmax nmax];
        else 
            % case: d = [0 -y z]
            % init guess is thrust -y z direction:
            % n = [nmax nmax nmax nmax].'
            % alpha = [-0.934 0 0.934 0]
            x0 = [0.934 0 -0.934 0 nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [0 y 0]
            % initial guess: max thrust y direction:
            x0 = [-pi/2 0 pi/2 0 nmax 133.3593 nmax 133.3593];
        elseif d(3) < 0
            % case: d = [0 y -z]
            % init guess is thrust y -z direction:
            x0 = [-pi/2 pi pi/2 pi nmax nmax nmax nmax];
        else 
            % case: d = [0 y z]
            % init guess is thrust y z direction:
            x0 = [-0.934 0 0.934 0 nmax nmax nmax nmax];
        end
    end
elseif d(1) < 0
    if d(2) == 0
        if d(3) == 0 
            % case: d = [-x 0 0]
            % initial guess: max thrust -x direction:
            x0 = [0 -pi/2 0 pi/2 133.3593 nmax 133.3593 nmax];
        elseif d(3) < 0 
            % case: d = [-x 0 -z]
            % initial guess: max thrust -x -z direction:
            x0 = [pi -pi/2 pi pi/2 nmax nmax nmax nmax];
        else
            % case: d = [-x 0 z]
            % init guess is max thrust -x z direction:
            x0 = [0 -0.934 0 0.934 nmax nmax nmax nmax];
        end
    elseif d(2) < 0
         if d(3) == 0 
            % case: d = [-x -y 0]
            % initial guess: max thrust in -x -y direction:
            x0 = [ acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2) acos(nhover^2/nmax^2) nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [-x -y -z]
            % init guess: max thrust in -x -y -z direction :
            x0 = [3*pi/4 -3*pi/4 -3*pi/4 3*pi/4 nmax nmax nmax nmax];
        else
            % case: d = [-x -y z]
            % init guess: max thrust in -x -y z direction :
            x0 = [pi/4 -pi/4 -pi/4 pi/4 nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [-x y 0]
            % initial guess: max thrust in x- y direction:
            x0 = [ -acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2), acos(nhover^2/nmax^2) acos(nhover^2/nmax^2) nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [-x y -z]
            % init guess: max thrust in -x y -z direction :
            x0 = [-3*pi/4 -3*pi/4 3*pi/4 3*pi/4 nmax nmax nmax nmax];
        else
            % case: d = [-x y z]
            % init guess: max thrust in -x y z direction :
            x0 = [-pi/4 -pi/4 pi/4 pi/4 nmax nmax nmax nmax];
        end
    end   
else
    if d(2) == 0
        if d(3) == 0 
            % case: d = [x 0 0]
            % initial guess: max thrust x direction:
            x0 = [0 0.934 0 -0.934 nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [x 0 -z]
            % initial guess: max thrust x -z direction:
            x0 = [pi pi/2 pi -pi/2 nmax nmax nmax nmax];
        else
            % case: d = [x 0 z]
            % init guess is max thrust x z direction:
            x0 = [0 0.934 0 -0.934 nmax nmax nmax nmax];
        end
    elseif d(2) < 0
         if d(3) == 0 
            % case: d = [x -y 0]
            % initial guess: max thrust in x -y direction:
           x0 = [acos(nhover^2/nmax^2) acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2) nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [x -y -z]
            % init guess: max thrust in x -y -z direction :
            x0 = [3*pi/4 3*pi/4 -3*pi/4 -3*pi/4 nmax nmax nmax nmax];
        else
            % case: d = [x -y z]
            % init guess: max thrust in x -y z direction :
            x0 = [pi/4 pi/4 -pi/4 -pi/4 nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [x y 0]
            % initial guess: max thrust in x y direction:
            x0 = [ -acos(nhover^2/nmax^2) acos(nhover^2/nmax^2), acos(nhover^2/nmax^2) -acos(nhover^2/nmax^2) nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [x y -z]
            % init guess: max thrust in x y -z direction :
            x0 = [-3*pi/4 3*pi/4 3*pi/4 -3*pi/4 nmax nmax nmax nmax];
        else
            % case: d = [x y z]
            % init guess: max thrust in x y z direction :
            x0 = [-pi/4 pi/4 pi/4 -pi/4 nmax nmax nmax nmax];
        end
    end
end

options = optimoptions('fmincon','Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
x = fmincon(@(x)fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconF(x, kf,d, m, g),  options);

% Substitution:
alpha = [x(1) x(2) x(3) x(4)];
n = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
k = 10.^Ndecimals;
n = round(k*n)/k;
alpha = round(k*alpha)/k;

end

function [c,ceq] = nonlinconF(x, kf,d, m, g)
%[c,ceq] = nonlincon(x, kf,d, m, g)
%NONLINCONF nonlinear constraint of the optimisation
% Hover condition
if d(3) < 0
    c = [];
else
    c(1) = -cos(x(1))*x(5)^2 - cos(x(2))*x(6)^2 - cos(x(3))*x(7)^2 - cos(x(4))*x(8)^2 +m*g/kf;
end
% Thrust parallel to d conditions: Ceq(x) = 0
ceq(1) = d(3)*(-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2)-d(2)*(cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2 -m*g/kf);
ceq(2) = -d(3)*(sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)+d(1)*(cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2 -m*g/kf);
ceq(3) = d(2)*(sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2)-d(1)*(-sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2);
end
