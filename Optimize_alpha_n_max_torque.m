function [alpha, n] = Optimize_alpha_n_max_torque(kf, km, nmax, g, d, m, L)
%[alpha, n] = Optimize_alpha_n_max_torque(kf, km, nmax, g, d, m, L)
%OPTIMIZE_ALPHA_N_MAX_TORQUE Optimization of alpha and n
%   maximize norm of the Torque M in an arbitrairy direction d:

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
            % initial guess: max torque -z direction:
            x0 = [1.164 1.164 1.164 1.164 nmax nmax nmax nmax];
            
        else 
            % case: d = [0 0 z]
            % initial guess: max torque z direction:
            x0 = [-1.164 -1.164 -1.164 -1.164 nmax nmax nmax nmax];
           
        end
    elseif d(2) < 0
        if d(3) == 0 
            % case: d = [0 -y 0]
            % initial guess: max torque y- direction:
            x0 = [0 0 pi 0 nmax 133.3593 nmax 133.3593];
        elseif d(3) < 0
            % case: d = [0 -y -z]
            % init guess is torque -y -z direction:
            x0 = [0.3 0.6 pi 0.6 nmax nmax nmax nmax];
        else 
            % case: d = [0 -y z]
            % init guess is thrust -y z direction:
            % n = [nmax nmax nmax nmax].'
            % alpha = [-0.934 0 0.934 0]
            x0 = [-0.3 -0.6 -pi -0.6 nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [0 y 0]
            % initial guess: max torque y direction:
            x0 = [pi 0 0 0 nmax 133.3593 nmax 133.3593];
        elseif d(3) < 0
            % case: d = [0 y -z]
            % init guess is torque y -z direction:
            x0 = [pi 0.6 0.3 0.6 nmax nmax nmax nmax];
        else 
            % case: d = [0 y z]
            % init guess is torque y z direction:
            x0 = [-pi -0.6 -0.3 -0.6 nmax nmax nmax nmax];
        end
    end
elseif d(1) < 0
    if d(2) == 0
        if d(3) == 0 
            % case: d = [-x 0 0]
            % initial guess: max torque -x direction:
            x0 = [0 pi 0 0 133.3593 nmax 133.3593 nmax];
        elseif d(3) < 0 
            % case: d = [-x 0 -z]
            % initial guess: max torque -x -z direction:
            x0 = [ 0.6 pi 0.6 0.3  nmax nmax nmax nmax];
        else
            % case: d = [-x 0 z]
            % init guess is max torque -x z direction:
            x0 = [ -0.6 -pi -0.6 -0.3  nmax nmax nmax nmax];
        end
    elseif d(2) < 0
         if d(3) == 0 
            % case: d = [-x -y 0]
            % initial guess: max torque in -x -y direction:
            x0 = [0 pi pi 0 nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [-x -y -z]
            % init guess: max thrust in -x -y -z direction :
            x0 = [0 pi pi 0 nmax nmax nmax nmax];
        else
            % case: d = [-x -y z]
            % init guess: max torque in -x -y z direction :
            x0 = [0 -pi -pi 0 nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [-x y 0]
            % initial guess: max torque in x- y direction:
            x0 = [pi pi 0 0 nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [-x y -z]
            % init guess: max torque in -x y -z direction :
            x0 = [pi pi 0 0 nmax nmax nmax nmax];
        else
            % case: d = [-x y z]
            % init guess: max torque in -x y z direction :
            x0 = [-pi -pi 0 0 nmax nmax nmax nmax];
        end
    end   
else
    if d(2) == 0
        if d(3) == 0 
            % case: d = [x 0 0]
            % initial guess: max torque x direction:
            x0 = [0 0 0 pi  133.3593 nmax 133.3593 nmax];
        elseif d(3) < 0 
            % case: d = [x 0 -z]
            % initial guess: max torque x -z direction:
            x0 = [0.6 0.3 0.6 pi nmax nmax nmax nmax];
        else
            % case: d = [x 0 z]
            % init guess is max torque x z direction:
            x0 = [-0.6 -0.3 -0.6 -pi nmax nmax nmax nmax];
        end
    elseif d(2) < 0
         if d(3) == 0 
            % case: d = [x -y 0]
            % initial guess: max torque in x -y direction:
            x0 = [0 0 pi pi nmax nmax nmax nmax];
         elseif d(3) < 0 
            % case: d = [x -y -z]
            % init guess: max torque in x -y -z direction :
            x0 = [0 0 pi pi nmax nmax nmax nmax];
        else
            % case: d = [x -y z]
            % init guess: max torque in x -y z direction :
            x0 = [0 0 -pi -pi nmax nmax nmax nmax];
        end
    else
        if d(3) == 0 
            % case: d = [x y 0]
            % initial guess: max torque in x y direction:
            x0 = [pi 0 0 pi nmax nmax nmax nmax];
        elseif d(3) < 0 
            % case: d = [x y -z]
            % init guess: max torque in x y -z direction :
            x0 = [pi 0 0 pi nmax nmax nmax nmax];
        else
            % case: d = [x y z]
            % init guess: max torque in x y z direction :
            x0 = [-pi 0 0 -pi nmax nmax nmax nmax];
        end
    end
end

options = optimoptions('fmincon','Algorithm','sqp');
options=optimoptions(options, 'MaxFunEvals',100000);
options=optimoptions(options,'MaxIter',100000);
% options = optimoptions(options, 'fmincon','Display','iter','Algorithm','sqp');
x = fmincon(@(x)fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlinconM(x, kf, km,d, m, g, L),  options);

% Substitution:
alpha = [x(1) x(2) x(3) x(4)];
n = [x(5) x(6) x(7) x(8)].';
Ndecimals = 6;
dec = 10.^Ndecimals;
n = round(dec*n)/dec;
alpha = round(dec*alpha)/dec;

end

function [c,ceq] = nonlinconM(x, kf,km, d, m, g, L)
%[c,ceq] = nonlincon(x, kf, nmax,d, m, g)
%NONLINCONM nonlinear constraint of the optimisation
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
ceq(6) = cos(x(1))*x(5)^2 + cos(x(2))*x(6)^2 + cos(x(3))*x(7)^2 + cos(x(4))*x(8)^2 -m*g/kf;
end

