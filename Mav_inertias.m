function [m, I] = Mav_inertias(n, L, theta, beta)
% computes masse of a n-rotor mav
% numerical values derived from voliro (with 6 rotors) 
mm2m = 10^(-6);
%% Mass of the drone body
mb = 2.166562421*n/6; % make the drone's body mass a fct of the number of propeller

%% Inertia of the drone body
Ibodyx = 7500.003560409*n/6*mm2m; % make the drone's body inertia a fct of the number of propeller
Ibodyy = 10938.848193227*n/6*mm2m; 
Ibodyz = 13694.808991418*n/6*mm2m;
Ibody = [Ibodyx, 0, 0; 0, Ibodyy, 0; 0, 0, Ibodyz];

%% Mass of an arm
mt = 0.07*L/0.2; % make mass of the tube (as fct of the length)
mp = 0.276873455-0.07; % mass of the propeller
ma = mt+ mp; % total mass of the arm

%% Inertia of an arm modeled as a tube of length L starting at the origin
r1 = 0.0045; % inner radius of the tube (measured on volro)
r2 = 0.0055; % outer radius of the tube (measured on volro)
Itube = [mt*(r1^2+r2^2)/2, 0 , 0; ...
        0, mt*(3*(r1^2+r2^2) + 4*L^2)/12, 0; ...
        0, 0, mt*(3*(r1^2+r2^2) + 4*L^2)/12]; % inertia of a tube

%% calculate the total inertia for all the arms
interval = 2*pi/n; % interval between arms in normal n-copter configuration
Iarms = zeros(3,3);
for i = 1:n
    Rb = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i));
    Ip = mp*(L^2*eye(3) - [L 0 0]*[L 0 0].'); % Model the propellers inertia as a point mass at length L of the center
    Iarms = Iarms + Rb*(Itube+Ip)*Rb.' ; % Add the inertias of the prop and the tube (in the body frame)
end

%% Total inertia and mass of the MAV
I = Ibody + Iarms;
m = ma + mb;
end