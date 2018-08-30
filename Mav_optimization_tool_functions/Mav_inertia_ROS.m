function [Irotorunit, mrotorunit, m, I] = Mav_inertia_ROS(n,beta, theta, L)
%%%% This file allows to calculate the inertia for the different design of
%%%% drones That are modeled on ROS

%% Parameters
mm2m = 10^(-6);

%% Mass of the drone body (based on values from voliro)
mb = 2.166562421*n/6;

%% Inertia of the drone body (based on values from voliro)
Ibx = 7500.003560409*n/6*mm2m;
Iby = 10938.848193227*n/6*mm2m; 
Ibz = 13694.808991418*n/6*mm2m;
Ibyz = -3.882565649*n/6*mm2m;
Ibxz =  24.695463703*n/6*mm2m; 
Ibxy = -34.208492595*n/6*mm2m;
Ib = [Ibx,  Ibxy, Ibxz; ...
      Ibxy, Iby,  Ibyz; ...
      Ibxz, Ibyz, Ibz];
  
%% Mass of an arm
mtspecifict = 0.1; % [kg/m]
ma = mtspecifict*L;
mp = (0.276873455-mtspecifict*0.3)/2;


%% Inertia tensor of the arm modeled as the inertia tensor of a tube   
r1 = 0.0045; % tube inner radius 
r2 = 0.0065; % tube outer radius 
Ia = [ma*(r1^2+r2^2)/2, 0 , 0; ...
          0, ma*(3*(r1^2+r2^2) + L^2)/12, 0; ...
          0, 0, ma*(3*(r1^2+r2^2) + L^2)/12];
r = [-L/2, 0, 0].'; % distance (from the center of the tube) at which we want the Inertia tensor
Ia = Ia + ma*(norm(r)^2*eye(3)-r*r.'); % parallel axis theorem (Steiner's rule)

%% Inertia of a propeller block modeled as a rectangle parallelepiped
w = 0.03; % widths
h = 0.065; % height
d = 0.03; % depth
Ip = [mp*(h^2+d^2)/12, 0 , 0; ...
        0, mp*((h^2+w^2))/12, 0; ...
        0, 0, mp*(3*(d^2+w^2))/12];
r = [0, 0, -h/2].'; % distance (from the center of the rect. parallel.) at which we want the Inertia tensor
Ip = Ip + mp*(norm(r)^2*eye(3)-r*r.'); % parallel axis theorem (Steiner's rule)

Irotorunit = Ip + Ia;
mrotorunit = ma + mp;

[m, I] = Mav_inertias(n, L, theta, beta);
end