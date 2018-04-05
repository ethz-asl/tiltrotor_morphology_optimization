%% Parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
m = Mb+4*Mp;% drone mass [kg]
g = 9.81;
 kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
% km = 2e-5;% Propeller drag coefficient
% L = 0.15; % Arm length [m]
nmax =150; % [roun/s]
nhover = sqrt((m*g/4)/kf); % [roun/s]
Ndecimals = 4;
dec = 10.^Ndecimals;
syms kf km L;
%% Simulation
syms teta1 teta2 teta3 teta4;
teta = [teta1 teta2 teta3 teta4];
syms beta1 beta2 beta3 beta4;
beta = [beta1 beta2 beta3 beta4];
syms alpha1 alpha2 alpha3 alpha4;
alpha = [alpha1 alpha2 alpha3 alpha4];
syms n1 n2 n3 n4;
n = [n1 n2 n3 n4];
 %teta = [0 0 0 0];
% beta = [0 0 0 0];
 %alpha = [0 0 0 0];
%% Rot mx
bRp1 = Rotz(teta(1))*Roty(beta(1))*Rotx(alpha(1));
bRp2 = Rotz(pi/2+teta(2))*Roty(beta(2))*Rotx(alpha(2));
bRp3 = Rotz(pi+teta(3))*Roty(beta(3))*Rotx(alpha(3));
bRp4 = Rotz(3*pi/2+teta(4))*Roty(beta(4))*Rotx(alpha(4));
Op1 = Rotz(teta(1))*Roty(beta(1))*[L 0 0].';
Op2 = Rotz(pi/2+teta(2))*Roty(beta(2))*[L 0 0].';
Op3 = Rotz(pi+teta(3))*Roty(beta(3))*[L 0 0].';
Op4 = Rotz(3*pi/2+teta(4))*Roty(beta(4))*[L 0 0].';

Tp1 = [0 0 kf*n1^2].';
Tp2 = [0 0 kf*n2^2].';
Tp3 = [0 0 kf*n3^2].';
Tp4 = [0 0 kf*n4^2].';

F  = bRp1*Tp1 + bRp2*Tp2 + bRp3*Tp3 + bRp4*Tp4;

Tauext1 = [0 0 km*n1^2].';
Tauext2 = [0 0 -km*n2^2].';
Tauext3 = [0 0 km*n3^2].';
Tauext4 = [0 0 -km*n4^2].';
% neglect rotor counter torque:
% Tauext1 = [0 0 0].';
% Tauext2 = [0 0 0].';
% Tauext3 = [0 0 0].';
% Tauext4 = [0 0 0].';

M = bRp1*Tauext1 + bRp2*Tauext2 + bRp3*Tauext3 + bRp4*Tauext4 + cross(Op1, bRp1*Tp1) + cross(Op2, bRp2*Tp2) + cross(Op3, bRp3*Tp3) + cross(Op4, bRp4*Tp4);
[row columns] = size(n);
A = [];
Astatic = [];
for ii = 1:3
    for jj = 1:columns
        Aij = F(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , {1,1,0}); Aij = subs(Aij,n , [0 0 0 0]);
        Ai2j = F(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)} , {1,1,pi/2}); Ai2j = subs(Ai2j,n , [0 0 0 0]);
        A = [A, Aij, Ai2j];
    end
     Astatic = [Astatic; A];
     A = [];
end
Astatic

A = [];
A_M_stat = [];
for ii = 1:3
    for jj = 1:columns
        Aij = M(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , {1,1,0});
        Aij = subs(Aij,n , [0 0 0 0]); Aij = subs(Aij,km , km/kf);
        Ai2j = M(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)} , {1,1,pi/2});
        Ai2j = subs(Ai2j,n , [0 0 0 0]); Ai2j = subs(Ai2j,km , km/kf);
        A = [A, Aij, Ai2j];
    end
     A_M_stat = [A_M_stat; A];
     A = [];
end
A_M_stat
% Astati
% Astatic = [sin(beta(1))*cos(teta(1)), sin(teta(1)), sin(beta(2))*cos(teta(2) + pi/2), sin(teta(2) + pi/2), -sin(beta(3))*cos(teta(3)), -sin(teta(3)), sin(beta(4))*cos(teta(4) + (3*pi)/2), sin(teta(4) + (3*pi)/2); ...
%            sin(beta(1))*sin(teta(1)), -cos(teta(1)), sin(beta(2))*sin(teta(2) + pi/2), -cos(teta(2) + pi/2), -sin(beta(3))*sin(teta(3)), cos(teta(3)), sin(beta(4))*sin(teta(4) + (3*pi)/2), -cos(teta(4) + (3*pi)/2); ...
%            cos(beta(1)), 0, cos(beta(2)), 0, cos(beta(3)), 0, cos(beta(4)),0];
%Fdec = [kf*cos(alpha(1))*n1^2; kf*sin(alpha(1))*n1^2; kf*cos(alpha(2))*n2^2; kf*sin(alpha(2))*n2^2; kf*cos(alpha(3))*n3^2; kf*sin(alpha(3))*n3^2; kf*cos(alpha(4))*n4^2; kf*sin(alpha(4))*n4^2];


% A_M_stat = [(km/kf*sin(beta(1))*cos(teta(1)) + L*sin(teta(1))),                      (km/kf*sin(teta(1)) - L*sin(beta(1))*cos(teta(1))), ...
%             (-km/kf*sin(beta(2))*cos(teta(2)+pi/2) + L*sin(teta(2)+pi/2)),           (-km/kf*sin(teta(2)+pi/2) - L*sin(beta(2))*cos(teta(2)+pi/2)), ...
%             (-km/kf*sin(beta(3))*cos(teta(3)) - L*sin(teta(3))),                     (-km/kf*sin(teta(3)) + L*sin(beta(3))*cos(teta(3))), ...
%             (-km/kf*sin(beta(4))*cos(teta(4)+(3*pi)/2) + L*sin(teta(4)+(3*pi)/2)),   (-km/kf*sin(teta(4)+(3*pi)/2) - L*sin(beta(4))*cos(teta(4)+(3*pi)/2)); ...
%             (km/kf*sin(beta(1))*sin(teta(1)) - L*cos(teta(1))),                      (-km/kf*cos(teta(1)) - L*sin(beta(1))*sin(teta(1))), ...
%             (-km/kf*sin(beta(2))*sin(teta(2)+pi/2) - L*cos(teta(2)+pi/2)),           (km/kf*cos(teta(2)+pi/2) - L*sin(beta(2))*sin(teta(2)+pi/2)), ...
%             (-km/kf*sin(beta(3))*sin(teta(3)) + L*cos(teta(3))),                     (km/kf*cos(teta(3)) + L*sin(beta(3))*sin(teta(3))), ...
%             (-km/kf*sin(beta(4))*sin(teta(4)+(3*pi)/2) - L*cos(teta(4)+(3*pi)/2)),   (km/kf*cos(teta(4)+(3*pi)/2) -L*sin(beta(4))*sin(teta(4)+(3*pi)/2)); ...
%             km/kf*cos(beta(1)), -L*cos(beta(1)), -km/kf*cos(beta(2)),
%             -L*cos(beta(2)), km/kf*cos(beta(3)),  -L*cos(beta(3)), -km/kf*cos(beta(4)), -L*cos(beta(4))];

% Astatic = [ sin(beta(1))*cos(teta(1)),  sin(teta(1)), sin(beta(2))*cos(teta(2) + pi/2),  sin(teta(2) + pi/2), …
% -sin(beta(3))*cos(teta(3)), -sin(teta(3)), sin(beta(4))*cos(teta(4) + (3*pi)/2),  sin(teta(4) + (3*pi)/2); …
% sin(beta(1))*sin(teta(1)), -cos(teta(1)), sin(beta(2))*sin(teta(2) + pi/2), -cos(teta(2) + pi/2), …
% -sin(beta(3))*sin(teta(3)),  cos(teta(3)), sin(beta(4))*sin(teta(4) + (3*pi)/2), -cos(teta(4) + (3*pi)/2); …
% cos(beta(1)),  0, cos(beta(2)), 0, cos(beta(3)), 0,  cos(beta(4)), 0];
% 
% A_M_stat = [ L*sin(teta(1))+(km*sin(beta(1))*cos(teta(1)))/kf, (km*sin(teta(1)))/kf-L*sin(beta(1))*cos(teta(1)),   …
% L*sin(teta(2) + pi/2)-(km*sin(beta(2))*cos(teta(2) + pi/2))/kf, -L*sin(beta(2))*cos(teta(2) + pi/2)-(km*sin(teta(2) + pi/2))/kf, …
% -L*sin(teta(3))-(km*sin(beta(3))*cos(teta(3)))/kf, L*sin(beta(3))*cos(teta(3))-(km*sin(teta(3)))/kf, …
% L*sin(teta(4) + (3*pi)/2)-(km*sin(beta(4))*cos(teta(4) + (3*pi)/2))/kf, -L*sin(beta(4))*cos(teta(4) + (3*pi)/2)-(km*sin(teta(4) + (3*pi)/2))/kf; …
%  (km*sin(beta(1))*sin(teta(1)))/kf-L*cos(teta(1)), -L*sin(beta(1))*sin(teta(1))-(km*cos(teta(1)))/kf, …
% -L*cos(teta(2) + pi/2)-(km*sin(beta(2))*sin(teta(2) + pi/2))/kf, (km*cos(teta(2) + pi/2))/kf-L*sin(beta(2))*sin(teta(2) + pi/2), …
% L*cos(teta(3))-(km*sin(beta(3))*sin(teta(3)))/kf, L*sin(beta(3))*sin(teta(3))+(km*cos(teta(3)))/kf, …
% -L*cos(teta(4) + (3*pi)/2)-(km*sin(beta(4))*sin(teta(4) + (3*pi)/2))/kf, (km*cos(teta(4) + (3*pi)/2))/kf-L*sin(beta(4))*sin(teta(4) + (3*pi)/2); ...
% (km*cos(beta(1)))/kf, - L*cos(beta(1)), -(km*cos(beta(2)))/kf, - L*cos(beta(2))*, …
% (km*cos(beta(3)))/kf, - L*cos(beta(3)),  -(km*cos(beta(4)))/kf, - L*cos(beta(4))];
