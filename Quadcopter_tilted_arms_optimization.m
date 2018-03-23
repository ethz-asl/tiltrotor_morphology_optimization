function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
g = 9.81;
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
nmax =150; % [roun/s]
nhover = sqrt(((Mb+4*Mp)*g/4)/kf); % [roun/s]
%% Inertia

%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb = rotz(roll0)*roty(pitch0)*rotz(yaw0);

%% Optimization of alpha and n 
d = [0 0 0].';
D = d;
% find optimal alpha and n for a max thrust in direction d
[alphastar, nstar] = Optimize_F_alpha_n(kf, nmax, nhover, g, d, m);
% calculate angular and linear acceleration with this alphastar and nstar
[pdotdot, wbdot] = quadRotorDynamic(kf, km, Ib, wRb, alphastar,nstar, L, g, m);
F = m*pdotdot; % Force produced by the MAV
Feff = 0;
Meff = 0;
% find optimal alpha and n for a max torque in direction d
[alphastar, nstar] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L);
% calculate angular and linear acceleration with this alphastar and nstar
[pdotdot, wbdot] = quadRotorDynamic(kf, km, Ib, wRb, alphastar,nstar, L, g, m);
M = Ib*wbdot; % Torque produced by the MAV

% Loop to compute the optimal Force in "any" directions:
for i = -1:1:1
    for j = -1:1:1
        for k = -1:1:1
            d = [i j k].';
            D = [D d];
            % find optimal alpha and n for a max torque in direction d
            [alphastar, nstar] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L);
            % calculate angular and linear acceleration with this alphastar and nstar
            [pdotdot, wbdot] = quadRotorDynamic(kf, km, Ib, wRb, alphastar,nstar, L, g, m);
            % Round the Force and Torque to the 4th decimal.
            Ndecimals = 4;
            k = 10.^Ndecimals;
            wbdot = round(k*wbdot)/k;
            M = [M Ib*wbdot];% Torque produced by the MAV
            Meff = [Meff L*kf*nstar.'*nstar];
            
            % find optimal alpha and n for a max thrust in direction d
            [alphastar, nstar] = Optimize_F_alpha_n(kf, nmax, nhover, g, d, m);
            % calculate angular and linear acceleration with this alphastar and nstar
            [pdotdot, wbdot] = quadRotorDynamic(kf, km, Ib, wRb, alphastar,nstar, L, g, m);
            pdotdot = round(k*pdotdot)/k;
            F = [F m*pdotdot];% Force produced by the MAV
            Feff = [Feff kf*nstar.'*nstar];
            
        end
    end
end
D_norm = vecnorm(D);
i0 = find(~D_norm);
D_unit = D./D_norm;
[C,ia,ic] = unique(D_unit.', 'stable', 'rows');
D = D(:,ia.');
F = F(:,ia.');
Feff = Feff(:,ia.');
Meff = Meff(:,ia.');
M = M(:,ia.');
i0 = find(~vecnorm(D));
D(:,i0) = [];
F(:,i0) = [];
Feff(:,i0) = [];
Meff(:,i0) = []; 
M(:,i0) = [];
%H = [D; F; M];
Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Fmax = max(vecnorm(F));
Fmin = min(vecnorm(F));
Mmax = max(vecnorm(M));
Mmin = min(vecnorm(M));


%% functions
function [alpha, n] = Optimize_F_alpha_n(kf, nmax, nhover, g, d, m)
% function [alpha, n] = Optimize_F_alpha_n(kf, nmax, nhover, g, d, m)

%% Optimization of alpha and n
% maximize norm of the Thrust F in an arbitrairy direction d:
% x = [alpha1 alpha2 alpha3 alpha4 n1 n2 n3 n4]
% => T = kf*[sin(x(2))*x(6)^2 - sin(x(4))*x(8)^2 ; -sin(x(1))*x(5)^2 + sin(x(3))*x(7)^2; ...
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
% function [c,ceq] = nonlincon(x, kf,d, m, g)
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

function [alpha, n] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L)
% function [alpha, n] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L)

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
k = 10.^Ndecimals;
n = round(k*n)/k;
alpha = round(k*alpha)/k;

end

function [c,ceq] = nonlinconM(x, kf,km, d, m, g, L)
% function [c,ceq] = nonlincon(x, kf, nmax,d, m, g)

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
end

