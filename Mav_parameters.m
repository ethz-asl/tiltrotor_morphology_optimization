function [g, dec, kf, km, alphamin, alphamax, wmin, wmax] = Mav_parameters()
%%%%%%%%%%%% Functoion containing all the parameters needed for the optimization %%%%%%%%%%%%
%% Parameters
g = 9.81;
Ndecimals = 5;
dec = 10.^Ndecimals; % Decimal to which rounnd the results
kf = 3.86e-4; % Propeller thrust coefficient % false[kg.m]
km = 1.5e-5;% Propeller drag coefficient
alphamin = -pi; % Minimum tilting angle
alphamax = pi; % Maximum tilting angle
wmin = 0; % minimum rotor speed allowed [round/s]
wmax =130; % maximum rotor speed allowed [round/s]
end