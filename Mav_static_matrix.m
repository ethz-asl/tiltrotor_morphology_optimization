function [A_F_static, A_M_static] = Mav_static_matrix(kf, km, L, beta, theta, n, dec)
%MAV_STATIC_MATRIX for a drone design with a 
% The static matrix are static allocation matrix that do not depend on the
% rotor orientation and speed.

% Vector containing all the decomposed vertical and horizontal forces:
% Fdec = [kf*cos(alpha(1))*w(1)^2; kf*sin(alpha(1))*w(1)^2;
%         kf*cos(alpha(2))*w(2)^2; kf*sin(alpha(2))*w(2)^2;
%         ...
%         kf*cos(alpha(n))*w(n)^2; kf*sin(alpha(n))*w(n)^2]; 

% Thrust in the body frame is: T = A_T_static*Fdec

% Torque in the body frame is: M = A_M_static*Fdec

%% Build Static matrix
syms km_sim;
alpha = sym('alpha',[1 n]);% Simulate alpha (propeller tilting angle)
w = sym('w',[1 n]);% Simulate w (propeller rotation speed [round/s])
interval = 2*pi/n; % interval between arms in normal n-copter configuration

% pre-allocation:
Taub = [0; 0; 0];
M = [0; 0; 0];
T = [0; 0; 0];
for i =1:n
    
    %% Find the propellers frame rotation matrix
    bRp(:,:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i))*Rotx(alpha(i));
    
    %% Find the propellers positions in the body frame
    Op(:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i))*[L 0 0].';
    
    %% Find forces applied by all propellers thrusts on the body
    Tp(:,i) = [0 0 w(i)^2].'; % force applied by every propeller in propeller frame
    T  = T + bRp(:,:,i)*Tp(:,i); % force applied by all the propellers in body frame
    
    %% Find torques applied by all propellers on the body (in propeller frame)
    if mod(i,2) == 0
        c = 1; % counter torque defined negative for clock wise rotating propellers (all odd propellers)
    else
        c = -1;
    end
    Tauext(:,i) = [0 0 c*km_sim*w(i)^2].'; % counter torque produced by every propeller in propeller frame
    M = M + bRp(:,:,i)*Tauext(:,i); % torque applied by every propeller on the body
    Taub = Taub + cross(Op(:,i),bRp(:,:,i)*Tp(:,i)); % total torque applied on the body (in body frame)
end
M = M + Taub;
%% Construct the force static matix
for ii = 1:3
    for jj = 1:n
        A_vertical_ij = T(ii);
        A_vertical_ij = subs(A_vertical_ij,{w(jj), alpha(jj)} , {1,0});
        A_vertical_ij = subs(A_vertical_ij,w , zeros(size(w))); 
        A_horizontal_ij = T(ii);
        A_horizontal_ij = subs(A_horizontal_ij,{w(jj), alpha(jj)}, {1,pi/2});
        A_horizontal_ij = subs(A_horizontal_ij,w , zeros(size(w)));
        A_F_static(ii,2*jj-1) = A_vertical_ij;
        A_F_static(ii,2*jj) = A_horizontal_ij;
    end
end
A_F_static = double(A_F_static);
A_F_static = round(A_F_static*dec)/dec;

%% Construct the torque static matix
for ii = 1:3
    for jj = 1:n
        A_vertical_ij = M(ii);
        A_vertical_ij = subs(A_vertical_ij,{w(jj), alpha(jj)} , {1,0});
        A_vertical_ij = subs(A_vertical_ij,w , zeros(size(w)));
        A_vertical_ij = subs(A_vertical_ij,km_sim , km/kf);
        A_horizontal_ij = M(ii);
        A_horizontal_ij = subs(A_horizontal_ij,{w(jj), alpha(jj)} , {1,pi/2});
        A_horizontal_ij = subs(A_horizontal_ij,w , zeros(size(w)));
        A_horizontal_ij = subs(A_horizontal_ij,km_sim , km/kf);
        A_M_static(ii,2*jj-1) = A_vertical_ij;
        A_M_static(ii,2*jj) = A_horizontal_ij;
    end
end
A_M_static = double(A_M_static);
A_M_static = round(A_M_static*dec)/dec;
end

