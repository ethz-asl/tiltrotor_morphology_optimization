function [g, dec, kf, km, step, max_iterations, optimize_alpha] = Mav_parameters()
%%%%%%%%%%%% Functoion containing all the parameters needed for the optimization %%%%%%%%%%%%
%% Parameters
g = 9.81;
Ndecimals = 5;
dec = 10.^Ndecimals; % 
kf = 3.86e-4; % Propeller thrust coefficient % false[kg.m]
km = 1.5e-5;% Propeller drag coefficient
%% Parameters for the optimization of tilting angles (alpha) and rotor speeds (w):
step = .1; % step to define the number of directions in which to compute forcetorque/hover eff
           % (0.5 -> 98 directions, 0.25 -> 578 directions, 0.1 -> 7490 directions)
max_iterations = 150; % Maximal number of times fmincom is iterated in one diection to find maximal force/maximal torque/ optimal hover mode
optimize_alpha = false; % If true performs an optimization on the tilting angles and rotor speeds to max the force/torque/hover eff in every direction
                       % if false uses the angles returned by the static matrix solution
                       
end
