function [F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot, optim)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust.m(beta ,theta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot)
%QUADCOPTER_MAX_TORQUE compute the maximal thrusts and torues for an
%arbitrary design of quadcopter 
%   Design defined by the design number ii, the arms angle (beta & theta)
%   and other parameters

%%%%%%%%%%%% Quadcopter with tilting rotor and tilted arms design optimization%%%%%%%%%%%%
%% Parameters
m = Mb+4*Mp;% drone mass [kg]
nhover = sqrt((m*g/4)/kf); % [roun/s]
Ndecimals = 4;
dec = 10.^Ndecimals;
%% init
roll0 = 0;
pitch0 = 0;
yaw0 = 0;
wRb = rotz(roll0)*roty(pitch0)*rotz(yaw0);
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

% Test
% Star = zeros(16,1);
% Fcheck= [0 0 0].';
% MM = zeros(3,1);

% Loop to compute the optimal Force in "any" directions (neglecting gravity):
for i = -1:step:1
    for j = -1:step:1
        for k = -1:step:1
            d = [i j k].';
            D = [D d];
            d0 = [0 0 0].';
            if isequal(d,d0) == 1
                d = [0;0;0];
            else
                d = d/vecnorm(d);
            end
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
            n0 = round(dec*n0)/dec;
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alpha0 = round(dec*alpha0)/dec;
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max thrust in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_thrust(kf, nmax, nhover, alpha0, n0, d, beta, theta);
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha0, beta, theta,n0, L, g, Mb, Mp, R);
            pdotdot = round(dec*pdotdot)/dec;
            F = [F m*pdotdot];% Force produced by the MAV
            Feff = [Feff kf*n0.'*n0];
            
%             % Double check            
%             Fdec_check = A_F_staticinv*(m*pdotdot);
%             nstar1 = [1/kf*sqrt(Fdec_check(1)^2 + Fdec_check(2)^2); 1/kf*sqrt(Fdec_check(3)^2 + Fdec_check(4)^2); 1/kf*sqrt(Fdec_check(5)^2 + Fdec_check(6)^2); 1/kf*sqrt(Fdec_check(7)^2 + Fdec_check(8)^2)];
%             nstar1 = sqrt(nstar1);
%             nstar1 = round(dec*nstar1)/dec;
%             alphastar1 = [atan2(Fdec_check(2),Fdec_check(1)) atan2(Fdec_check(4),Fdec_check(3)) atan2(Fdec_check(6),Fdec_check(5)) atan2(Fdec_check(8),Fdec_check(7))];
%             alphastar1 = round(dec*alphastar1)/dec;
%             if norm(alphastar1) < 4*pi
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar1, beta, theta,nstar1, L, g, Mb, Mp, R);
%                 pdotdot = round(dec*pdotdot)/dec;
%             else
%                 pdotdot =[0 0 0].';
%             end
%             Fcheck = [Fcheck m*pdotdot];
%             star = [nstar;alphastar.';nstar1;alphastar1.'];
%             Star = [Star star];
             
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
            n0 = round(dec*n0)/dec;
            alpha0 = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alpha0 = round(dec*alpha0)/dec;
            if ~isequal(d, [0 0 0].')
                if optim
                    % Perform the optimization and find the max torque in direction d 
                    [alphastar, nstar] = Quadcopter_tilted_arms_max_torque(kf, km, L, nmax, nhover, alpha0, n0, d, beta, theta);                 
                else
                    alphastar = alpha0;
                    nstar = n0;
                end
            else
                alphastar = [0 0 0 0];
                nstar = [0 0 0 0].';
            end
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha0, beta, theta,n0, L, g, Mb, Mp, R);
            pdotdot = round(dec*pdotdot)/dec;
            
            M = [M Ib*wbdot];% Force produced by the MAV
            Meff = [Meff L*kf*n0.'*n0];
            
%             % Test
%             star = [nstar;alphastar.'];
%             Star = [Star star];
%             MM = [MM Mdes];
        end
    end
end
D_norm = vecnorm(D);
i0 = find(~D_norm);
D_norm(:,i0) = [];
D(:,i0) = [];
F(:,i0) = [];
Feff(:,i0) = [];
M(:,i0) = [];
Meff(:,i0) = [];


% % Test
% Fcheck(:,i0) = [];
% Star(:,i0) = [];
% MM(:,i0) = [];

D_unit = D./D_norm;
[C,ia,ic] = unique(D_unit.', 'stable', 'rows');
D = D(:,ia.');
F = F(:,ia.');
Feff = Feff(:,ia.');
Meff = Meff(:,ia.');
M = M(:,ia.');

% % Test
% Fcheck = Fcheck(:,ia.');
% Star = Star(:,ia.');
% MM = MM(:,ia.');


Feff = 100*vecnorm(F)./Feff;
Meff = 100*vecnorm(M)./Meff;
Fnorm = vecnorm(F);
Fmax = max(Fnorm);
Fmin = min(Fnorm);
Mnorm = vecnorm(M);
Mmax = max(Mnorm);
Mmin = min(Mnorm);

% % Test
% MM = [M;MM;Star];
if plot
figure(ii); 
subplot(2,2,1);
colormap(flipud(jet));
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

    % Generate a sphere consisting of 20by 20 faces of radius
    [x,y,z]=sphere;
    % use surf function to plot
    hSurface=surf(3*x,3*y,3*z);
    hold on;
    set(hSurface,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none');
end
daspect([1 1 1]);
title(['Max torque space for design ' num2str(ii) ''])
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

subplot(2,2,2);
colormap(flipud(jet));
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
    % Generate a sphere consisting of 20by 20 faces of radius Fmin
    [x,y,z]=sphere;
    % use surf function to plot
    hSurface=surf(x,y,z);
    hold on
    set(hSurface,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none')
    %axis([-20 20 -20 20 -20 20]);
end
daspect([1 1 1]);
title(['Max torque space for design ' num2str(ii) ''])
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight



subplot(2,2,3);
% Plot drone design
% First Simple design representation in grey (beta = theta = 0)

% Generate a sphere
[x,y,z]=sphere;
% use surf function to plot
R= R/3;
r = 2*R/5;
if ~isequal(beta,[0 0 0 0])
    propelerSphere1=surf(r*x+L,r*y,r*z); hold on;% centered at Op1
    set(propelerSphere1,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    propelerSphere2=surf(r*x+0,r*y+L,r*z+0); % centered at Op2
    set(propelerSphere2,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    propelerSphere3=surf(r*x-L,r*y,r*z); % centered at Op3
    set(propelerSphere3,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    propelerSphere4=surf(r*x,r*y-L,r*z); % centered at Op4
    set(propelerSphere4,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 L], [0 0], [0 0],'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
    plot3([0 0], [0 L], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
    plot3([0 -L], [0 0], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
    plot3([0 0], [0 -L], [0 0], 'Color',[0 0 0]+0.05*k, 'LineWidth', 40*R)
end

%The actual design with tilted arms
centerSphere=surf(R*x,R*y,R*z);% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
propelerSphere1=surf(r*x+Op1(1),r*y+Op1(2),r*z+Op1(3)); % centered at Op1
% text(Op1(1)+0.02*sign(Op1(1)), Op1(2)+0.02*sign(Op1(2)), Op1(3)+0.01*sign(Op1(3)),'Op1')
set(propelerSphere1,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere2=surf(r*x+Op2(1),r*y+Op2(2),r*z+Op2(3)); % centered at Op2
% text(Op2(1)+0.02*sign(Op2(1)), Op2(2)+0.02*sign(Op2(2)), Op2(3)+0.01*sign(Op2(3)),'Op2')
set(propelerSphere2,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere3=surf(r*x+Op3(1),r*y+Op3(2),r*z+Op3(3)); % centered at Op3
% text(Op3(1)+0.02*sign(Op3(1)), Op3(2)+0.02*sign(Op3(2)), Op3(3)+0.01*sign(Op3(3)),'Op3')
set(propelerSphere3,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere4=surf(r*x+Op4(1),r*y+Op4(2),r*z+Op4(3)); % centered at Op4
% text(Op4(1)+0.02*sign(Op4(1)), Op4(2)+0.02*sign(Op4(2)), Op4(3)+0.01*sign(Op4(3)),'Op4')
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

x0=10;
y0=10;
width=1020;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height]);
daspect([1 1 1]);
title(['Illustration of design ' num2str(ii) ])
str = (['Design ' num2str(ii) ': \beta = [' num2str(rad2deg(beta(1))) ', ' ...
         num2str(rad2deg(beta(2))) ', ' num2str(rad2deg(beta(3))) ...
         ', ' num2str(rad2deg(beta(4))) '], \theta = [' ...
         num2str(round(rad2deg(theta(1)))) ', ' num2str(round(rad2deg(theta(2)))) ', ' ...
         num2str(round(rad2deg(theta(3)))) ', ' num2str(round(rad2deg(theta(4)))) ...
         ']']);
dim = [ .4 .7 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlabel('X') 
ylabel('Y')
zlabel('Z')
camlight
end
end
