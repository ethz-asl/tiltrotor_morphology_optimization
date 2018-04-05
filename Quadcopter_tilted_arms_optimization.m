clear all;
close all;
%% Design
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
g = 9.81;
kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
km = 2e-5;% Propeller drag coefficient
L = 0.15; % Arm length [m]
nmax =150; % [roun/s]
beta = [0 0 0 0];
teta = [0 0 0 0];
ii = 1;
step = 1;
plot = false;

%% Find design properties:

[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopte_max_torque_thrust(pi/10*[0 0 0 0] ,teta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, true);

% for i = 0:pi/10:pi
%     [D,Fcheck, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopte_max_torque_thrust([i -i -i i] ,teta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, true);
%     A1 = [ii, rad2deg(i), rad2deg(-i), rad2deg(-i), rad2deg(i), Fmin, Fmax, Mmin, Mmax];                 
%     formatSpec = 'Design %3.0f (beta = [%3.0f %3.0f %3.0f %3.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f\n';
% 	fprintf(formatSpec,A1)
%     ii = ii+1;
% end

% beta_step = 1;
% Data = [];
% for i = -1:beta_step:1
%     for j = -1:beta_step:1
%         for k = -1:beta_step:1
%             for l = -1:beta_step:1
%                 if ii ==41 || ii== 25 || ii == 57
%                     plot= true;
%                 else 
%                     plot = false;
%                 end
%                 beta = pi/6*[i j k l];
%                 [D, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopte_max_torque_thrust(beta ,teta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot);
%                 A1 = [ii, rad2deg(beta(1)), rad2deg(beta(2)), rad2deg(beta(3)), rad2deg(beta(4)), Fmin, Fmax, Mmin, Mmax];
%                 formatSpec = 'Design %3.0f (beta = [%3.0f %3.0f %3.0f %3.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f\n';
%                 fprintf(formatSpec,A1)
%                 ii = ii + 1;
%                 Data = [Data, [Fmin;Fmax;Mmin;Mmax;beta(1);beta(2);beta(3);beta(4)]];
%             end
%         end
%     end
% end
