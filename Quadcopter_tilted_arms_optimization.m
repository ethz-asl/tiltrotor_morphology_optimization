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
ii = 1;
step = 0.5;
plot = true;
optimize_theta = false;
optimize_max_thrust_torque = true;

%% Find design properties:

%[D,F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot, optimize_max_thrust_torque);

beta_step = 1;
Data = [];
for i = -1:beta_step:1
    for j = -1:beta_step:1
        for k = -1:beta_step:1
            for l = -1:beta_step:1
                if ii ==41 || ii== 25 || ii == 57 || ii ==75 || ii== 73 || ii == 52 || ii == 30 || ii ==9 || ii== 7
                    plot= true;
                else 
                    plot = false;
                end
                beta = pi/6*[i j k l];
                theta = Quadcopter_tilted_arms_find_theta(beta, optimize_theta);
                [F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot, optimize_max_thrust_torque);
                A1 = [ii, rad2deg(beta(1)), rad2deg(beta(2)), rad2deg(beta(3)), rad2deg(beta(4)),rad2deg(theta(1)), rad2deg(theta(2)), rad2deg(theta(3)), rad2deg(theta(4)), Fmin, Fmax, Mmin, Mmax];
                formatSpec = 'Design %3.0f (beta = [%3.0f %3.0f %3.0f %3.0f], theta = [%3.0f %3.0f %3.0f %3.0f]) -> Fmin = %2.2f, Fmax = %2.2f, Mmin = %2.2f Mmax = %2.2f\n';
                fprintf(formatSpec,A1)
                ii = ii + 1;
                Data = [Data, [Fmin;Fmax;Mmin;Mmax;beta(1);beta(2);beta(3);beta(4)]];
            end
        end
    end
end
