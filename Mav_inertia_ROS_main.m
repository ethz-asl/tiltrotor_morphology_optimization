%%%% This file allows to calculate the inertia for the different design of
%%%% drones that are modeled on ROS and use these inertias in ROS
close all;
clear all;
[file_path] = fileparts(mfilename('fullpath'));
addpath(file_path);
addpath([file_path '/Mav_optimization_tool_functions/']);

%% Enter the parameters of your drone (number of prop, angles beta, angles theta, arm length)
n=8;
if mod(n,2) ~= 0
	beta =  0.615479709*ones(1,n);
else
    for i = 1:n
        if mod(i,2) ~= 0
            beta(i) =  0.615479709;
        else
            beta(i) =  -0.615479709;
        end
    end
end
theta = zeros(1,n);

L = 0.5;

%% Will return you the inertia and mass of the drone 
%% and the inertia and mass of a propeller unit (with the arm)
[Irotorunit, mrotorunit, m, I] = Mav_inertia_ROS(n,beta, theta, L)
