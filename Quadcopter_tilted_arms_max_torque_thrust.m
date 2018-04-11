function [D, Heff, Hmin, Hmax, F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot, optim)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust.m(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot)
%QUADCOPTER_MAX_TORQUE compute the maximal thrusts and torues for an
%arbitrary design of quadcopter 
%   Design defined by the design number ii, the arms angle (beta & theta)
%   and other parameters

%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
m = Mb+4*Mp;% drone mass [kg]
Ndecimals = 4;
dec = 10.^Ndecimals;
%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb0 = rotz(rad2deg(roll0))*roty(rad2deg(pitch0))*rotz(rad2deg(yaw0));
%% Static matrix definition
% Fdec = [kf*cos(alpha(1))*n1^2; kf*sin(alpha(1))*n1^2; kf*cos(alpha(2))*n2^2; 
%         kf*sin(alpha(2))*n2^2; kf*cos(alpha(3))*n3^2; kf*sin(alpha(3))*n3^2; 
%         kf*cos(alpha(4))*n4^2; kf*sin(alpha(4))*n4^2];

% F = m*p'' = A_F_static*Fdec
% M = Ib*wb' = A_M_static*Fdec

% These static matrix are found using the file: Quadcopter_tilted_arms_Find_Static_Matrix.m
A_F_static = [sin(beta(1))*cos(theta(1)),  sin(theta(1)), sin(beta(2))*cos(theta(2) + pi/2),  sin(theta(2) + pi/2), ...
           -sin(beta(3))*cos(theta(3)), -sin(theta(3)), sin(beta(4))*cos(theta(4) + (3*pi)/2),  sin(theta(4) + (3*pi)/2); ...
           sin(beta(1))*sin(theta(1)), -cos(theta(1)), sin(beta(2))*sin(theta(2) + pi/2), -cos(theta(2) + pi/2), ...
           -sin(beta(3))*sin(theta(3)),  cos(theta(3)), sin(beta(4))*sin(theta(4) + (3*pi)/2), -cos(theta(4) + (3*pi)/2); ...
           cos(beta(1)),  0, cos(beta(2)), 0, cos(beta(3)), 0,  cos(beta(4)), 0];

A_M_static = [L*sin(theta(1))+(km*sin(beta(1))*cos(theta(1)))/kf, (km*sin(theta(1)))/kf-L*sin(beta(1))*cos(theta(1)),   ...
            L*sin(theta(2) + pi/2)-(km*sin(beta(2))*cos(theta(2) + pi/2))/kf, -L*sin(beta(2))*cos(theta(2) + pi/2)-(km*sin(theta(2) + pi/2))/kf, ...
            -L*sin(theta(3))-(km*sin(beta(3))*cos(theta(3)))/kf, L*sin(beta(3))*cos(theta(3))-(km*sin(theta(3)))/kf, ...
            L*sin(theta(4) + (3*pi)/2)-(km*sin(beta(4))*cos(theta(4) + (3*pi)/2))/kf, -L*sin(beta(4))*cos(theta(4) + (3*pi)/2)-(km*sin(theta(4) + (3*pi)/2))/kf; ...
            (km*sin(beta(1))*sin(theta(1)))/kf-L*cos(theta(1)), -L*sin(beta(1))*sin(theta(1))-(km*cos(theta(1)))/kf, ...
            -L*cos(theta(2) + pi/2)-(km*sin(beta(2))*sin(theta(2) + pi/2))/kf, (km*cos(theta(2) + pi/2))/kf-L*sin(beta(2))*sin(theta(2) + pi/2), ...
            L*cos(theta(3))-(km*sin(beta(3))*sin(theta(3)))/kf, L*sin(beta(3))*sin(theta(3))+(km*cos(theta(3)))/kf, ...
            -L*cos(theta(4) + (3*pi)/2)-(km*sin(beta(4))*sin(theta(4) + (3*pi)/2))/kf, (km*cos(theta(4) + (3*pi)/2))/kf-L*sin(beta(4))*sin(theta(4) + (3*pi)/2); ...
            (km*cos(beta(1)))/kf, - L*cos(beta(1)), -(km*cos(beta(2)))/kf, - L*cos(beta(2)), ...
            (km*cos(beta(3)))/kf, - L*cos(beta(3)),  -(km*cos(beta(4)))/kf, - L*cos(beta(4))];

% => Fdec = inv(A_F_static)*Fdes
A_F_staticinv = pinv(A_F_static);

% => Fdec = inv(A_M_static)*Mdes
A_M_staticinv = pinv(A_M_static);

%% Optimization of alpha and n 
% initialization:
d = [0 0 0].';
D = d;
F = [0 0 0].';
Feff = 0;
M = zeros(3,1);
Meff = 0;
Heff = 0;

% Test
errF = 0;
errM = 0;
errH = 0;
% Star = [];
% FH = zeros(3,1);
% StarOpt = zeros(10,1);
% Fcheck= [0 0 0].';
% MM = zeros(3,1);

%% Loop to compute the optimal Force in "any" directions (neglecting gravity):
for i = -1:step:1
    for j = -1:step:1
        for k = -1:step:1
            d = [i j k].';
            D = [D d];
            d0 = [0 0 0].';
            if isequal(d,d0)
                d = [0;0;0];
            else
                d = d./vecnorm(d);
            end
            %% find max force in direction d
            % find initial alpha and n for the optimisation to find the max thrust in direction d
            Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers
            Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
            % Inverse substitution : 
            %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
            %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
            n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
            n0 = n0/vecnorm(n0);
            n0 = nmax^2*n0/max(n0); % 0 <= nstar <= nmax
            n0(n0<0) = 0; 
            n0 = sqrt(n0);
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
             % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta,n0, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            F0 = m*pdotdot;
            FN0 = norm(F0);
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max thrust in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, alpha0, n0, d, beta, theta);
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta,nstar, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Fstar = m*pdotdot;
            FNstar = norm(Fstar);
            if FNstar < FN0
                F = [F F0];% Force produced by the MAV
                Feff = [Feff kf*n0.'*n0];
            else
                F = [F Fstar];% Force produced by the MAV
                Feff = [Feff kf*nstar.'*nstar];
            end
%             %test
%             if isequal(alphastar, alpha0) && isequal(nstar, n0)
%             	errF = errF+1;
%             end

%             % Double check            
%             Fdec_check = A_F_staticinv*(m*pdotdot);
%             nstar1 = [1/kf*sqrt(Fdec_check(1)^2 + Fdec_check(2)^2); 1/kf*sqrt(Fdec_check(3)^2 + Fdec_check(4)^2); 1/kf*sqrt(Fdec_check(5)^2 + Fdec_check(6)^2); 1/kf*sqrt(Fdec_check(7)^2 + Fdec_check(8)^2)];
%             nstar1 = sqrt(nstar1);
%             nstar1 = round(dec*nstar1)/dec;
%             alphastar1 = [atan2(Fdec_check(2),Fdec_check(1)) atan2(Fdec_check(4),Fdec_check(3)) atan2(Fdec_check(6),Fdec_check(5)) atan2(Fdec_check(8),Fdec_check(7))];
%             alphastar1 = round(dec*alphastar1)/dec;
%             if norm(alphastar1) < 4*pi
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar1, beta, theta,nstar1, L, g, Mb, Mp, R, false);
%                 pdotdot = round(dec*pdotdot)/dec;
%             else
%                 pdotdot =[0 0 0].';
%             end
%             Fcheck = [Fcheck m*pdotdot];
%             star = [nstar;alphastar.';nstar1;alphastar1.'];
%             Star = [Star star];
            
            %% find max torque in direction d 
            % find initial alpha and n for the optimisation to find the max torque in direction d
            Mdes = d*(4*L*nmax^2*kf); % set desired torque to be equal to the maximal torque of the four propellers
            Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
            % Inverse substitution : 
            %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
            %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
            n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
            n0 = n0/vecnorm(n0);
            n0 = nmax^2*n0/max(n0); % 0 <= nstar <= nmax
            n0(n0<0) = 0; 
            n0 = sqrt(n0);
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alpha0, beta, theta, n0, L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            M0 = Ib*wbdot;
            MN0 = norm(M0);
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max torque in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, alpha0, n0, d, beta, theta);                 
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta, nstar, L, g, Mb, Mp, R, false);
            wbdot = round(dec*wbdot)/dec;
            Mstar = Ib*wbdot;
            MNstar = norm(Mstar);
            if MNstar < MN0
                M = [M M0];% Force produced by the MAV
                Meff = [Meff L*kf*n0.'*n0];
            else
                M = [M Mstar];% Force produced by the MAV
                Meff = [Meff L*kf*nstar.'*nstar];
            end
%             % Test
%             if isequal(alphastar, alpha0) && isequal(nstar, n0)
%             	errM = errM+1;
%             end
%             star = [nstar;alphastar.'];
%             Star = [Star, star];

            %% find hover efficiency in direction d 
            if ~isequal(d, [0 0 0].')
                roll = 0;
                pitch = acos(d(3));
                yaw = atan2(d(2),d(1));
                wRb = rotz(rad2deg(roll))*roty(rad2deg(pitch))*rotz(rad2deg(yaw));
                % find initial alpha and n for the optimisation to find the max torque in direction d
                Fdes = m*g*d;
                Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
                % Inverse substitution : 
                %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
                %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
                n0 = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
                n0 = sqrt(n0);
                alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
                %Test
%                 star = [round(n0); round(rad2deg(alpha0)).'; 0; round(100*m*g/(kf*norm(n0)^2))];
%                 Star = [Star star];
                if optim
                    % Perform the optimization and find the max torque in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_opt_hover(kf, Fdes, nmax, alpha0, n0, beta, theta);
                    nstar = round(dec*nstar)/dec;
                    n0 = round(dec*n0)/dec;
                    Hstar = m*g/(kf*norm(nstar)^2);
                    H0 = m*g/(kf*norm(n0)^2);
                    if Hstar < H0
                        Heff = [Heff H0]; % Hover efficiency in direction d
                    else
                        Heff = [Heff Hstar]; % Hover efficiency in direction d
                    end
                else
                    alphastar = alpha0;
                    alphastar = round(dec*alphastar)/dec;
                    nstar = n0;
                    nstar = round(dec*nstar)/dec;
                    Hstar = m*g/(kf*norm(nstar)^2);
                    Heff = [Heff Hstar]; % Hover efficiency in direction d
                end
                %Test
%                 if isequal(alphastar, alpha0) && isequal(nstar, n0)
%                     errH = errH+1;
%                 end
%                 alphastar = round(dec*alphastar)/dec;
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb0, alphastar, beta, theta,nstar, L, g, Mb, Mp, R, false);
%                 pdotdot = round(dec*pdotdot)/dec;
%                 FH = [FH (m*pdotdot-m*g*d)];                
%                 staropt = [round(nstar); round(rad2deg(alphastar)).'; 0; round(100*(m*g/(kf*norm(nstar)^2)))];
%                 StarOpt = [StarOpt staropt];
                
            else
                Heff = [Heff 0];
            end
        end
    end
end
i0 = find(~vecnorm(D));
D(:,i0) = [];
F(:,i0) = [];
Feff(:,i0) = [];
M(:,i0) = [];
Meff(:,i0) = [];
Heff(:,i0) = [];

% % Test
% Fcheck(:,i0) = [];
% Star(:,i0) = [];
% FH(:,i0) = [];
% StarOpt(:,i0) = [];
% MM(:,i0) = [];
% errF
% errM
% errH

D_unit = D./vecnorm(D);
[D_unit,ia,ic] = unique(D_unit.', 'stable', 'rows');
D = D(:,ia.');
F = F(:,ia.');
Feff = Feff(:,ia.');
Meff = Meff(:,ia.');
Heff= Heff(:,ia.');
M = M(:,ia.');

% % Test
% Fcheck = Fcheck(:,ia.');
% FH = FH(:,ia.');
% Star = Star(:,ia.');
% StarOpt = StarOpt(:,ia.');
% MM = MM(:,ia.');


Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Heff = 100*Heff;
Heff = round(Heff*dec/100)/(dec/100);
Fnorm = vecnorm(F);
Fmax = max(Fnorm);
Fmin = min(Fnorm);
Mnorm = vecnorm(M);
Mmax = max(Mnorm);
Mmin = min(Mnorm);
Hmax = max(Heff);
Hmin = min(Heff);

% % Test
% MM = [M;MM;Star];
if plot
figure(ii); 
%% Plot Force space
subplot(2,2,1);
colormap(flipud(jet));
%Feff = round(Feff);
scatter3(F(1,:), F(2,:), F(3,:),  100 ,Feff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of the Thrust (%)';
caxis([min(Feff) max(Feff)]);
% Link the neighbours with one another (draws the edges of the force space polyedron)
if step >= 0.5
    neighbour1 = [D(1,:); D(2,:)];
    neighbour2 = [D(2,:); D(3,:)];
    neighbour3 = [D(1,:); D(3,:)];
    [row column] = size(D);
    for jj = 1:column
        for kk = 1:column
            if jj~=kk
                if ~isequal(D(:,jj),-D(:,kk))
                    if neighbour1(:,jj) == neighbour1(:,kk)
                        if abs(D(3,jj)-D(3,kk))<2*step
                            plot3([F(1,jj) F(1,kk)], [F(2,jj) F(2,kk)], [F(3,jj) F(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour2(:,jj) == neighbour2(:,kk)
                        if abs(D(1,jj)-D(1,kk))<2*step
                            plot3([F(1,jj) F(1,kk)], [F(2,jj) F(2,kk)], [F(3,jj) F(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour3(:,jj) == neighbour3(:,kk)
                        if abs(D(2,jj)-D(2,kk))<2*step
                            plot3([F(1,jj) F(1,kk)], [F(2,jj) F(2,kk)], [F(3,jj) F(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end    
                end
            end
        end
    end
end
daspect([1 1 1]);
title('Reachable force space')
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% Plot torque space
subplot(2,2,2);
colormap(flipud(jet));
%Meff = round(Meff);
scatter3(M(1,:), M(2,:), M(3,:),  100 ,Meff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of the Torque (%)';
caxis([min(Meff) max(Meff)])
% Link the neighbours with one another (draws the edges of the torque space polyedron)
if step >= 0.5
    for jj = 1:column
        for kk = 1:column
            if jj~=kk
                if ~isequal(D(:,jj),-D(:,kk))
                    if neighbour1(:,jj) == neighbour1(:,kk)
                        if abs(D(3,jj)-D(3,kk))<2*step
                            plot3([M(1,jj) M(1,kk)], [M(2,jj) M(2,kk)], [M(3,jj) M(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour2(:,jj) == neighbour2(:,kk)
                        if abs(D(1,jj)-D(1,kk))<2*step
                            plot3([M(1,jj) M(1,kk)], [M(2,jj) M(2,kk)], [M(3,jj) M(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour3(:,jj) == neighbour3(:,kk)
                        if abs(D(2,jj)-D(2,kk))<2*step
                            plot3([M(1,jj) M(1,kk)], [M(2,jj) M(2,kk)], [M(3,jj) M(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end    
                end
            end
        end
    end
end
daspect([1 1 1]);
title('Reachable torque space')
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% Plot Hover efficincy
subplot(2,2,3);
colors = flipud(jet);
colormap(colors);
%Heff = round(Heff);
D_unit = D_unit.';
scatter3(D_unit(1,:), D_unit(2,:), D_unit(3,:),  100 ,Heff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of Hover (%)';
caxis([Hmin Hmax])
% Link the neighbours with one another (draws the edges of the torque space polyedron)
if step >= 0.5
    for jj = 1:column
        for kk = 1:column
            if jj~=kk
                if ~isequal(D(:,jj),-D(:,kk))
                    if neighbour1(:,jj) == neighbour1(:,kk)
                        if abs(D(3,jj)-D(3,kk))<2*step
                            plot3([D_unit(1,jj) D_unit(1,kk)], [D_unit(2,jj) D_unit(2,kk)], [D_unit(3,jj) D_unit(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour2(:,jj) == neighbour2(:,kk)
                        if abs(D(1,jj)-D(1,kk))<2*step
                            plot3([D_unit(1,jj) D_unit(1,kk)], [D_unit(2,jj) D_unit(2,kk)], [D_unit(3,jj) D_unit(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end
                    if neighbour3(:,jj) == neighbour3(:,kk)
                        if abs(D(2,jj)-D(2,kk))<2*step
                            plot3([D_unit(1,jj) D_unit(1,kk)], [D_unit(2,jj) D_unit(2,kk)], [D_unit(3,jj) D_unit(3,kk)],'Color',[0.5 0.5 0.5], 'LineWidth', 1);
                        end
                    end    
                end
            end
        end
    end
end
daspect([1 1 1]);
title('Hover efficiency')
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% Plot drone representation
subplot(2,2,4);
% Generate a sphere
[x,y,z]=sphere;
% use surf function to plot
R= R/3;
r = 2*R/5;
centerSphere=surf(R*x,R*y,R*z);% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
propelerSphere1=surf(r*x+Op1(1),r*y+Op1(2),r*z+Op1(3)); % centered at Op1
if ~isequal(Op1, [L 0 0].')
    % plot the original quadcopter design 
    propelerSphere1=surf(r*x+L,r*y,r*z); hold on;% centered at Op1
    set(propelerSphere1,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 L], [0 0], [0 0],'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
end
set(propelerSphere1,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere2=surf(r*x+Op2(1),r*y+Op2(2),r*z+Op2(3)); % centered at Op2
if ~isequal(Op1, [0 L 0].')
    % plot the original quadcopter design 
    propelerSphere2=surf(r*x+0,r*y+L,r*z+0); % centered at Op2
    set(propelerSphere2,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 0], [0 L], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
end
set(propelerSphere2,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere3=surf(r*x+Op3(1),r*y+Op3(2),r*z+Op3(3)); % centered at Op3
if ~isequal(Op1, [-L 0 0].')
    % plot the original quadcopter design 
    propelerSphere3=surf(r*x-L,r*y,r*z); % centered at Op3
    set(propelerSphere3,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 -L], [0 0], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
end
set(propelerSphere3,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere4=surf(r*x+Op4(1),r*y+Op4(2),r*z+Op4(3)); % centered at Op4
if ~isequal(Op1, [0 -L 0].')
    % plot the original quadcopter design 
    propelerSphere4=surf(r*x,r*y-L,r*z); % centered at Op4
    set(propelerSphere4,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 0], [0 -L], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
end
set(propelerSphere4,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
plot3([0 Op1(1)], [0 Op1(2)], [0 Op1(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op2(1)], [0 Op2(2)], [0 Op2(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op3(1)], [0 Op3(2)], [0 Op3(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op4(1)], [0 Op4(2)], [0 Op4(3)], 'c', 'LineWidth', 400*R)
quiver3(0, 0, 0, 0, 0, -0.1, 'b')

% plot angle theta
if theta(1) ~= 0
    Op10 = roty(rad2deg(beta(1)))*[L 0 0].';
    Op10norm = norm(Op10);
    Op10(3) = 0;
    Op10 = cos(beta(1))*Op10*Op10norm/norm(Op10);
    Op1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*[L 0 0].';
    Op1norm = norm(Op1);
    Op1(3) = 0;
    Op1 = cos(beta(1))*Op1*Op1norm/norm(Op1);
    theta01 = [[0; 0; 0], Op10, Op1 ];
    fill3(theta01(1,:),theta01(2,:),theta01(3,:),'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op1+Op10)/2;
    text(postxt(1), postxt(2), postxt(3), '\theta_{1}')
end
if theta(2) ~= 0
    Op20 = rotz(rad2deg(pi/2))*roty(rad2deg(beta(2)))*[L 0 0].';
    Op20norm = norm(Op20);
    Op20(3) = 0;
    Op20 = cos(beta(2))*Op20*Op20norm/norm(Op20);
    Op2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*[L 0 0].';
    Op2norm = norm(Op2);
    Op2(3) = 0;
    Op2 = cos(beta(2))*Op2*Op2norm/norm(Op2);
    theta02 = [[0; 0; 0], Op20, Op2 ];
    fill3(theta02(1,:),theta02(2,:),theta02(3,:),'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op2+Op20)/2;
    text(postxt(1), postxt(2), postxt(3), '\theta_{2}')
end
if theta(3) ~= 0
    Op30 = rotz(rad2deg(pi))*roty(rad2deg(beta(3)))*[L 0 0].';
    Op30norm = norm(Op30);
    Op30(3) = 0;
    Op30 = cos(beta(3))*Op30*Op30norm/norm(Op30);
    Op3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*[L 0 0].';
    Op3norm = norm(Op3);
    Op3(3) = 0;
    Op3 = cos(beta(3))*Op3*Op3norm/norm(Op3);
    theta03 = [[0; 0; 0], Op30, Op3 ];
    fill3(theta03(1,:),theta03(2,:),theta03(3,:),'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op3+Op30)/2;
    text(postxt(1), postxt(2), postxt(3), '\theta_{3}')
end
if theta(3) ~= 0
    Op40 = rotz(rad2deg(3*pi/2))*roty(rad2deg(beta(4)))*[L 0 0].';
    Op40norm = norm(Op40);
    Op40(3) = 0;
    Op40 = cos(beta(4))*Op40*Op40norm/norm(Op40);
    Op4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*[L 0 0].';
    Op4(3) = 0;
    Op4norm = norm(Op4);
    Op4(3) = 0;
    Op4 = Op4*Op4norm/norm(Op4);
    theta04 = [[0; 0; 0], Op40, Op4 ];
    fill3(theta04(1,:),theta04(2,:),theta04(3,:),'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op4+Op40)/2;
    text(postxt(1), postxt(2), postxt(3), '\theta_{4}')
end


% plot angle beta
if beta(1) ~= 0
    Op10 = cos(beta(1))*rotz(rad2deg(theta(1)))*[L 0 0].';
    Op1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*[L 0 0].';
    beta01 = [[0; 0; 0], Op10, Op1 ];
    fill3(beta01(1,:),beta01(2,:),beta01(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op1+Op10)/2;
    text(postxt(1), postxt(2), postxt(3), '\beta_{1}') 
end
if beta(2) ~= 0
    Op20 = cos(beta(2))*rotz(rad2deg(pi/2+theta(2)))*[L 0 0].';
    Op2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*[L 0 0].';
    beta02 = [[0; 0; 0], Op20, Op2 ];
    fill3(beta02(1,:),beta02(2,:),beta02(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op2+Op20)/2;
    text(postxt(1), postxt(2), postxt(3), '\beta_{2}') 
end
if beta(3) ~= 0
    Op30 = cos(beta(3))*rotz(rad2deg(pi+theta(3)))*[L 0 0].';
    Op3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*[L 0 0].';
    beta03 = [[0; 0; 0], Op30, Op3 ];
    fill3(beta03(1,:),beta03(2,:),beta03(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op3+Op30)/2;
    text(postxt(1), postxt(2), postxt(3), '\beta_{3}') 
end
if beta(3) ~= 0
    Op40 = cos(beta(4))*rotz(rad2deg(3*pi/2+theta(4)))*[L 0 0].';
    Op4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*[L 0 0].';
    beta04 = [[0; 0; 0], Op40, Op4 ];
    fill3(beta04(1,:),beta04(2,:),beta04(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op4+Op40)/2;
    text(postxt(1), postxt(2), postxt(3), '\beta_{4}') 
end
daspect([1 1 1]);
title('Illustration of the design')
xlabel('X') 
ylabel('Y')
zlabel('Z')
camlight

%% General plot options and annotations
x0=10;
y0=10;
width=1020;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height]);
str = (['Design ' num2str(ii) ': \beta = [' num2str(rad2deg(beta(1))) ', ' ...
         num2str(rad2deg(beta(2))) ', ' num2str(rad2deg(beta(3))) ...
         ', ' num2str(rad2deg(beta(4))) '], \theta = [' ...
         num2str(round(rad2deg(theta(1)))) ', ' num2str(round(rad2deg(theta(2)))) ', ' ...
         num2str(round(rad2deg(theta(3)))) ', ' num2str(round(rad2deg(theta(4)))) ...
         '], ' 'F_{min} = ' num2str(Fmin) ', F_{max} = ' num2str(Fmax) ...
         ', M_{min} = ' num2str(Mmin) ', M_{max} = ' num2str(Mmax) ...
         ', H_{min} = ' num2str(Hmin) ', H_{max} = ' num2str(Hmax)]);
dim = [ .2 .7 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
end
end