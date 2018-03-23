%%%%%%%%%%%% Quadcopter with tilting rotor design optimization%%%%%%%%%%%%
clear all;
close all;
%% Parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
m = Mb + 4*Mp; % total mass [kg]
g = 9.81;
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
nmax =150; % [roun/s]
nhover = sqrt((m*g/4)/kf); % [roun/s]
%% Inertia
Op1 = [L 0 0].';
Op2 = [0 L 0].';
Op3 = [-L 0 0].';
Op4 = [0 -L 0].';
Icom = 2/5*Mb*R*R*eye(3);
% Inertia tensor of a sphere with rotors represented as point masses
Ib = Mp*(norm(Op1)^2*eye(3) - Op1*Op1.' + norm(Op2)^2*eye(3) - Op2*Op2.' + ...
     norm(Op3)^2*eye(3) - Op3*Op3.' + norm(Op4)^2*eye(3) - Op4*Op4.') +Icom;
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
[alphastar2, nstar2] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L);
% calculate angular and linear acceleration with this alphastar and nstar
[pdotdot2, wbdot2] = quadRotorDynamic(kf, km, Ib, wRb, alphastar2,nstar2, L, g, m);
M = Ib*wbdot; % Torque produced by the MAV
FM = m*pdotdot2;
% Loop to compute the optimal Force in "any" directions:
for i = -1:1:1
    for j = -1:1:1
        for k = -1:1:1
            d = [i j k].';
            D = [D d];
            % find optimal alpha and n for a max torque in direction d
            [alphastar2, nstar2] = Optimize_M_alpha_n(kf, km, nmax, g, d, m, L);
            % calculate angular and linear acceleration with this alphastar and nstar
            [pdotdot2, wbdot2] = quadRotorDynamic(kf, km, Ib, wRb, alphastar2,nstar2, L, g, m);
            % Round the Force and Torque to the 4th decimal.
            Ndecimals = 4;
            k = 10.^Ndecimals;
            wbdot2 = round(k*wbdot2)/k;
            M = [M Ib*wbdot2];% Torque produced by the MAV
            pdotdot2 = round(k*pdotdot2)/k;
            FM = [FM m*pdotdot2];
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
FM = FM(:,ia.');
Feff = Feff(:,ia.');
Meff = Meff(:,ia.');
M = M(:,ia.');
i0 = find(~vecnorm(D));
D(:,i0) = [];
F(:,i0) = [];
FM(:,i0) = [];
Feff(:,i0) = [];
Meff(:,i0) = []; 
M(:,i0) = [];
MF = [M; FM];
Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Fmax = max(vecnorm(F));
Fmin = min(vecnorm(F));
Mmax = max(vecnorm(M));
Mmin = min(vecnorm(M));

F1 = [F(:,14), F(:,17),  F(:,26), F(:,23)];% (0 0 1), (0 1 1),  (1 1 1), (1 0 1)
F2 = [F(:,14), F(:,23), F(:,20), F(:,12)];% (0 0 1), (1 0 1), (1 -1 1), (0 -1 1)
F3 = [F(:,14), F(:,6), F(:,3), F(:,12)];% (0 0 1), (-1 0 1), (-1 -1 1), (0 -1 1)
F4 = [F(:,14), F(:,6), F(:,9), F(:,17)];% (0 0 1), (-1 0 1), (-1 1 1), (0 1 1)
F5 = [F(:,22), F(:,23),  F(:,26), F(:,25)];% (1 0 0), (1 0 1),  (1 1 1), (1 1 0)
F6 = [F(:,22), F(:,23), F(:,20), F(:,19)];% (1 0 0), (1 0 1), (1 -1 1), (1-1 0)
F7 = [F(:,22), F(:,21), F(:,18), F(:,19)];% (1 0 0), (1 0 -1), (1 -1 -1), (1 -1 0)
F8 = [F(:,22), F(:,21), F(:,24), F(:,25)];% (1 0 0), (1 0 -1), (1 1 -1), (1 1 0)
F9 = [F(:,16), F(:,17),  F(:,26), F(:,25)];% (0 1 0), (0 1 1),  (1 1 1), (-1 1 0)
F10 = [F(:,16), F(:,17), F(:,9), F(:,8)];% (0 1 0), (0 1 1), (-1 1 1), (1 1 0)
F11 = [F(:,16), F(:,15), F(:,24), F(:,25)];% (0 1 0), (0 1 -1), (1 1 -1), (1 1 0)
F12 = [F(:,16), F(:,15), F(:,7), F(:,8)];% (0 1 0), (0 1 -1), (-1 1 -1), (-1 1 0)
F13 = [F(:,13), F(:,15),  F(:,24), F(:,21)];% (0 0 -1), (0 1 -1),  (1 1 -1), (1 0 -1)
F14 = [F(:,13), F(:,21), F(:,18), F(:,10)];% (0 0 -1), (1 0 -1), (1 -1 -1), (0 -1 -1)
F15 = [F(:,13), F(:,4), F(:,1), F(:,10)];% (0 0 -1), (-1 0 -1), (-1 -1 -1), (0 -1 -1)
F16 = [F(:,13), F(:,4), F(:,7), F(:,15)];% (0 0 -1), (-1 0 -1), (-1 1 -1), (0 1 -1)
F17 = [F(:,5), F(:,6),  F(:,9), F(:,8)];% (-1 0 0), (-1 0 1),  (-1 1 1), (-1 1 0)
F18 = [F(:,5), F(:,6), F(:,3), F(:,2)];% (-1 0 0), (-1 0 1), (-1 -1 1), (-1-1 0)
F19 = [F(:,5), F(:,4), F(:,1), F(:,2)];% (-1 0 0), (-1 0 -1), (-1 -1 -1), (-1 -1 0)
F20 = [F(:,5), F(:,4), F(:,7), F(:,8)];% (-1 0 0), (-1 0 -1), (-1 1 -1), (-1 1 0)
F21 = [F(:,11), F(:,12),  F(:,3), F(:,2)];% (0 -1 0), (0 -1 1),  (-1 -1 1), (-1 -1 0)
F22 = [F(:,11), F(:,12), F(:,20), F(:,19)];% (0 -1 0), (0 -1 1), (1 -1 1), (1 -1 0)
F23 = [F(:,11), F(:,10), F(:,18), F(:,19)];% (0 -1 0), (0 -1 -1), (1 -1 -1), (1 -1 0)
F24 = [F(:,11), F(:,10), F(:,1), F(:,2)];% (0 -1 0), (0 -1 -1), (-1 -1 -1), (-1 -1 0)

% z = 0 plan
Fz01 = [[0; 0; 0], F(:,11),  F(:,19), F(:,22)];% [0; 0; 0], (0 -1 0),  (1 -1 0), (1 0 0)
Fz02 = [[0; 0; 0], F(:,11),  F(:,2), F(:,5)];% [0; 0; 0], (0 -1 0),  (-1 -1 0), (-1 0 0)
Fz03 = [[0; 0; 0], F(:,16),  F(:,8), F(:,5)];% [0; 0; 0], (0 1 0),  (-1 1 0), (-1 0 0)
Fz04 = [[0; 0; 0], F(:,16),  F(:,25), F(:,22)];% [0; 0; 0], (0 1 0),  (1 1 0), (1 0 0)

figure(1); 
colormap(flipud(jet(20)));    
scatter3(F(1,:), F(2,:), F(3,:),  100 ,Feff, 'filled'); hold on;
c = colorbar;
c.Label.String = 'Efficiency of the Thrust (%)';

fill3(F1(1,:),F1(2,:),F1(3,:),'b', 'FaceAlpha', 0.2); hold on;
fill3(F2(1,:),F2(2,:),F2(3,:),'b', 'FaceAlpha', 0.2);
fill3(F3(1,:),F3(2,:),F3(3,:),'b', 'FaceAlpha', 0.2);
fill3(F4(1,:),F4(2,:),F4(3,:),'b', 'FaceAlpha', 0.2);
fill3(F5(1,:),F5(2,:),F5(3,:),'b', 'FaceAlpha', 0.2);
fill3(F6(1,:),F6(2,:),F6(3,:),'b', 'FaceAlpha', 0.2);
fill3(F7(1,:),F7(2,:),F7(3,:),'b', 'FaceAlpha', 0.2);
fill3(F8(1,:),F8(2,:),F8(3,:),'b', 'FaceAlpha', 0.2);
fill3(F9(1,:),F9(2,:),F9(3,:),'b', 'FaceAlpha', 0.2);
fill3(F10(1,:),F10(2,:),F10(3,:),'b', 'FaceAlpha', 0.2);
fill3(F11(1,:),F11(2,:),F11(3,:),'b', 'FaceAlpha', 0.2);
fill3(F12(1,:),F12(2,:),F12(3,:),'b', 'FaceAlpha', 0.2);
fill3(F13(1,:),F13(2,:),F13(3,:),'b', 'FaceAlpha', 0.2);
fill3(F14(1,:),F14(2,:),F14(3,:),'b', 'FaceAlpha', 0.2);
fill3(F15(1,:),F15(2,:),F15(3,:),'b', 'FaceAlpha', 0.2);
fill3(F16(1,:),F16(2,:),F16(3,:),'b', 'FaceAlpha', 0.2);
fill3(F17(1,:),F17(2,:),F17(3,:),'b', 'FaceAlpha', 0.2);
fill3(F18(1,:),F18(2,:),F18(3,:),'b', 'FaceAlpha', 0.2);
fill3(F19(1,:),F19(2,:),F19(3,:),'b', 'FaceAlpha', 0.2);
fill3(F20(1,:),F20(2,:),F20(3,:),'b', 'FaceAlpha', 0.2);
fill3(F21(1,:),F21(2,:),F21(3,:),'b', 'FaceAlpha', 0.2);
fill3(F22(1,:),F22(2,:),F22(3,:),'b', 'FaceAlpha', 0.2);
fill3(F23(1,:),F23(2,:),F23(3,:),'b', 'FaceAlpha', 0.2);
fill3(F24(1,:),F24(2,:),F24(3,:),'b', 'FaceAlpha', 0.2);
fill3(Fz01(1,:),Fz01(2,:),Fz01(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Fz02(1,:),Fz02(2,:),Fz02(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Fz03(1,:),Fz03(2,:),Fz03(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Fz04(1,:),Fz04(2,:),Fz04(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
% Generate a sphere consisting of 20by 20 faces of radius Fmin
[x,y,z]=sphere;
% use surf function to plot
hSurface=surf(3*x,3*y,3*z);
hold on;
set(hSurface,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none')
%axis([-20 20 -20 20 -20 20]);
daspect([1 1 1]);
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

M1 = [M(:,14), M(:,17),  M(:,26), M(:,23)];% (0 0 1), (0 1 1),  (1 1 1), (1 0 1)
M2 = [M(:,14), M(:,23), M(:,20), M(:,12)];% (0 0 1), (1 0 1), (1 -1 1), (0 -1 1)
M3 = [M(:,14), M(:,6), M(:,3), M(:,12)];% (0 0 1), (-1 0 1), (-1 -1 1), (0 -1 1)
M4 = [M(:,14), M(:,6), M(:,9), M(:,17)];% (0 0 1), (-1 0 1), (-1 1 1), (0 1 1)
M5 = [M(:,22), M(:,23),  M(:,26), M(:,25)];% (1 0 0), (1 0 1),  (1 1 1), (1 1 0)
M6 = [M(:,22), M(:,23), M(:,20), M(:,19)];% (1 0 0), (1 0 1), (1 -1 1), (1-1 0)
M7 = [M(:,22), M(:,21), M(:,18), M(:,19)];% (1 0 0), (1 0 -1), (1 -1 -1), (1 -1 0)
M8 = [M(:,22), M(:,21), M(:,24), M(:,25)];% (1 0 0), (1 0 -1), (1 1 -1), (1 1 0)
M9 = [M(:,16), M(:,17),  M(:,26), M(:,25)];% (0 1 0), (0 1 1),  (1 1 1), (-1 1 0)
M10 = [M(:,16), M(:,17), M(:,9), M(:,8)];% (0 1 0), (0 1 1), (-1 1 1), (1 1 0)
M11 = [M(:,16), M(:,15), M(:,24), M(:,25)];% (0 1 0), (0 1 -1), (1 1 -1), (1 1 0)
M12 = [M(:,16), M(:,15), M(:,7), M(:,8)];% (0 1 0), (0 1 -1), (-1 1 -1), (-1 1 0)
M13 = [M(:,13), M(:,15),  M(:,24), M(:,21)];% (0 0 -1), (0 1 -1),  (1 1 -1), (1 0 -1)
M14 = [M(:,13), M(:,21), M(:,18), M(:,10)];% (0 0 -1), (1 0 -1), (1 -1 -1), (0 -1 -1)
M15 = [M(:,13), M(:,4), M(:,1), M(:,10)];% (0 0 -1), (-1 0 -1), (-1 -1 -1), (0 -1 -1)
M16 = [M(:,13), M(:,4), M(:,7), M(:,15)];% (0 0 -1), (-1 0 -1), (-1 1 -1), (0 1 -1)
M17 = [M(:,5), M(:,6),  M(:,9), M(:,8)];% (-1 0 0), (-1 0 1),  (-1 1 1), (-1 1 0)
M18 = [M(:,5), M(:,6), M(:,3), M(:,2)];% (-1 0 0), (-1 0 1), (-1 -1 1), (-1-1 0)
M19 = [M(:,5), M(:,4), M(:,1), M(:,2)];% (-1 0 0), (-1 0 -1), (-1 -1 -1), (-1 -1 0)
M20 = [M(:,5), M(:,4), M(:,7), M(:,8)];% (-1 0 0), (-1 0 -1), (-1 1 -1), (-1 1 0)
M21 = [M(:,11), M(:,12),  M(:,3), M(:,2)];% (0 -1 0), (0 -1 1),  (-1 -1 1), (-1 -1 0)
M22 = [M(:,11), M(:,12), M(:,20), M(:,19)];% (0 -1 0), (0 -1 1), (1 -1 1), (1 -1 0)
M23 = [M(:,11), M(:,10), M(:,18), M(:,19)];% (0 -1 0), (0 -1 -1), (1 -1 -1), (1 -1 0)
M24 = [M(:,11), M(:,10), M(:,1), M(:,2)];% (0 -1 0), (0 -1 -1), (-1 -1 -1), (-1 -1 0)
% Z = 0 plan
Mz01 = [[0; 0; 0], M(:,11),  M(:,19), M(:,22)];% [0; 0; 0], (0 -1 0),  (1 -1 0), (1 0 0)
Mz02 = [[0; 0; 0], M(:,11),  M(:,2), M(:,5)];% [0; 0; 0], (0 -1 0),  (-1 -1 0), (-1 0 0)
Mz03 = [[0; 0; 0], M(:,16),  M(:,8), M(:,5)];% [0; 0; 0], (0 1 0),  (-1 1 0), (-1 0 0)
Mz04 = [[0; 0; 0], M(:,16),  M(:,25), M(:,22)];% [0; 0; 0], (0 1 0),  (1 1 0), (1 0 0)

figure(2);
colormap(flipud(jet(20)));    
scatter3(M(1,:), M(2,:), M(3,:),  100 ,Meff, 'filled'); hold on;
c = colorbar;
c.Label.String = 'Efficiency of the Torque (%)';

fill3(M1(1,:),M1(2,:),M1(3,:),'b', 'FaceAlpha', 0.4); hold on;
fill3(M2(1,:),M2(2,:),M2(3,:),'b', 'FaceAlpha', 0.4);
fill3(M3(1,:),M3(2,:),M3(3,:),'b', 'FaceAlpha', 0.4);
fill3(M4(1,:),M4(2,:),M4(3,:),'b', 'FaceAlpha', 0.4);
fill3(M5(1,:),M5(2,:),M5(3,:),'b', 'FaceAlpha', 0.4);
fill3(M6(1,:),M6(2,:),M6(3,:),'b', 'FaceAlpha', 0.4);
fill3(M7(1,:),M7(2,:),M7(3,:),'b', 'FaceAlpha', 0.4);
fill3(M8(1,:),M8(2,:),M8(3,:),'b', 'FaceAlpha', 0.4);
fill3(M9(1,:),M9(2,:),M9(3,:),'b', 'FaceAlpha', 0.4);
fill3(M10(1,:),M10(2,:),M10(3,:),'b', 'FaceAlpha', 0.4);
fill3(M11(1,:),M11(2,:),M11(3,:),'b', 'FaceAlpha', 0.4);
fill3(M12(1,:),M12(2,:),M12(3,:),'b', 'FaceAlpha', 0.4);
fill3(M13(1,:),M13(2,:),M13(3,:),'b', 'FaceAlpha', 0.4);
fill3(M14(1,:),M14(2,:),M14(3,:),'b', 'FaceAlpha', 0.4);
fill3(M15(1,:),M15(2,:),M15(3,:),'b', 'FaceAlpha', 0.4);
fill3(M16(1,:),M16(2,:),M16(3,:),'b', 'FaceAlpha', 0.4);
fill3(M17(1,:),M17(2,:),M17(3,:),'b', 'FaceAlpha', 0.4);
fill3(M18(1,:),M18(2,:),M18(3,:),'b', 'FaceAlpha', 0.4);
fill3(M19(1,:),M19(2,:),M19(3,:),'b', 'FaceAlpha', 0.4);
fill3(M20(1,:),M20(2,:),M20(3,:),'b', 'FaceAlpha', 0.4);
fill3(M21(1,:),M21(2,:),M21(3,:),'b', 'FaceAlpha', 0.4);
fill3(M22(1,:),M22(2,:),M22(3,:),'b', 'FaceAlpha', 0.4);
fill3(M23(1,:),M23(2,:),M23(3,:),'b', 'FaceAlpha', 0.4);
fill3(M24(1,:),M24(2,:),M24(3,:),'b', 'FaceAlpha', 0.4);
fill3(Mz01(1,:),Mz01(2,:),Mz01(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Mz02(1,:),Mz02(2,:),Mz02(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Mz03(1,:),Mz03(2,:),Mz03(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
fill3(Mz04(1,:),Mz04(2,:),Mz04(3,:),'k', 'FaceAlpha', 0.3, 'EdgeColor','none');

% Generate a sphere consisting of 20by 20 faces of radius Fmin
[x,y,z]=sphere;
% use surf function to plot
hSurface=surf(x,y,z);
hold on
set(hSurface,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none')
%axis([-20 20 -20 20 -20 20]);
daspect([1 1 1]);
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight


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