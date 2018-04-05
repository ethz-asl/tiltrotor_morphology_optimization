function [F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust.m(beta ,teta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot)
%[F,Fmin, Fmax, Feff, M, Mmin, Mmax, Meff] = Quadcopter_tilted_arms_max_torque_thrust.m(beta ,teta, ii, L, Mb, Mp, R, kf, km, nmax, g, step, plot)
%Quadcopter_tilted_arms_max_torque_thrust.m compute the maximal thrusts and torues for a
%design of quadcopter. design 
%   Quadcopter with tilting rotor and tilted arms d

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
%% Substitution
% Fdec = [kf*cos(alpha(1))*n1^2; kf*sin(alpha(1))*n1^2; kf*cos(alpha(2))*n2^2; 
%         kf*sin(alpha(2))*n2^2; kf*cos(alpha(3))*n3^2; kf*sin(alpha(3))*n3^2; 
%         kf*cos(alpha(4))*n4^2; kf*sin(alpha(4))*n4^2];

% F = m*p'' = A_F_static*Fdec
% M = Ib*wb' = A_M_static*Fdec

A_F_static = [sin(beta(1))*cos(teta(1)),  sin(teta(1)), sin(beta(2))*cos(teta(2) + pi/2),  sin(teta(2) + pi/2), ...
           -sin(beta(3))*cos(teta(3)), -sin(teta(3)), sin(beta(4))*cos(teta(4) + (3*pi)/2),  sin(teta(4) + (3*pi)/2); ...
           sin(beta(1))*sin(teta(1)), -cos(teta(1)), sin(beta(2))*sin(teta(2) + pi/2), -cos(teta(2) + pi/2), ...
           -sin(beta(3))*sin(teta(3)),  cos(teta(3)), sin(beta(4))*sin(teta(4) + (3*pi)/2), -cos(teta(4) + (3*pi)/2); ...
           cos(beta(1)),  0, cos(beta(2)), 0, cos(beta(3)), 0,  cos(beta(4)), 0];

A_M_static = [L*sin(teta(1))+(km*sin(beta(1))*cos(teta(1)))/kf, (km*sin(teta(1)))/kf-L*sin(beta(1))*cos(teta(1)),   ...
            L*sin(teta(2) + pi/2)-(km*sin(beta(2))*cos(teta(2) + pi/2))/kf, -L*sin(beta(2))*cos(teta(2) + pi/2)-(km*sin(teta(2) + pi/2))/kf, ...
            -L*sin(teta(3))-(km*sin(beta(3))*cos(teta(3)))/kf, L*sin(beta(3))*cos(teta(3))-(km*sin(teta(3)))/kf, ...
            L*sin(teta(4) + (3*pi)/2)-(km*sin(beta(4))*cos(teta(4) + (3*pi)/2))/kf, -L*sin(beta(4))*cos(teta(4) + (3*pi)/2)-(km*sin(teta(4) + (3*pi)/2))/kf; ...
            (km*sin(beta(1))*sin(teta(1)))/kf-L*cos(teta(1)), -L*sin(beta(1))*sin(teta(1))-(km*cos(teta(1)))/kf, ...
            -L*cos(teta(2) + pi/2)-(km*sin(beta(2))*sin(teta(2) + pi/2))/kf, (km*cos(teta(2) + pi/2))/kf-L*sin(beta(2))*sin(teta(2) + pi/2), ...
            L*cos(teta(3))-(km*sin(beta(3))*sin(teta(3)))/kf, L*sin(beta(3))*sin(teta(3))+(km*cos(teta(3)))/kf, ...
            -L*cos(teta(4) + (3*pi)/2)-(km*sin(beta(4))*sin(teta(4) + (3*pi)/2))/kf, (km*cos(teta(4) + (3*pi)/2))/kf-L*sin(beta(4))*sin(teta(4) + (3*pi)/2); ...
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
            % find optimal alpha and n for a max thrust in direction d
            Fdes = d.*(4*nmax^2*kf); % set desired force to be equal to the maximal thrust of the four propellers
            Fdec = A_F_staticinv*(Fdes); % Fdec = inv(Astatic)*Fdes
            % Inverse substitution : 
            %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
            %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
            nstar = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
            nstar = nstar/vecnorm(nstar);
            nstar = nmax^2*nstar/max(nstar); % 0 <= nstar <= nmax
            nstar(nstar<0) = 0; 
            nstar = sqrt(nstar);
            nstar = round(dec*nstar)/dec;
            alphastar = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alphastar = round(dec*alphastar)/dec;
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar, beta, teta,nstar, L, g, Mb, Mp, R);
            pdotdot = round(dec*pdotdot)/dec;
            F = [F m*pdotdot];% Force produced by the MAV
            Feff = [Feff kf*nstar.'*nstar];
            
%             % Double check            
%             Fdec_check = A_F_staticinv*(m*pdotdot);
%             nstar1 = [1/kf*sqrt(Fdec_check(1)^2 + Fdec_check(2)^2); 1/kf*sqrt(Fdec_check(3)^2 + Fdec_check(4)^2); 1/kf*sqrt(Fdec_check(5)^2 + Fdec_check(6)^2); 1/kf*sqrt(Fdec_check(7)^2 + Fdec_check(8)^2)];
%             nstar1 = sqrt(nstar1);
%             nstar1 = round(dec*nstar1)/dec;
%             alphastar1 = [atan2(Fdec_check(2),Fdec_check(1)) atan2(Fdec_check(4),Fdec_check(3)) atan2(Fdec_check(6),Fdec_check(5)) atan2(Fdec_check(8),Fdec_check(7))];
%             alphastar1 = round(dec*alphastar1)/dec;
%             if norm(alphastar1) < 4*pi
%                 [m, Ib,pdotdot, wbdot] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar1, beta, teta,nstar1, L, g, Mb, Mp, R);
%                 pdotdot = round(dec*pdotdot)/dec;
%             else
%                 pdotdot =[0 0 0].';
%             end
%             Fcheck = [Fcheck m*pdotdot];
%             star = [nstar;alphastar.';nstar1;alphastar1.'];
%             Star = [Star star];
             
            % find optimal alpha and n for a max torque in direction d
            Mdes = d*(4*L*nmax^2*kf); % set desired torque to be equal to the maximal torque of the four propellers
            Fdec = A_M_staticinv*(Mdes); % Fdec = inv(Astatic)*Fdes
            % Inverse substitution : 
            %                       ni² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
            %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
            nstar = [1/kf*sqrt(Fdec(1)^2 + Fdec(2)^2); 1/kf*sqrt(Fdec(3)^2 + Fdec(4)^2); 1/kf*sqrt(Fdec(5)^2 + Fdec(6)^2); 1/kf*sqrt(Fdec(7)^2 + Fdec(8)^2)];
            nstar = nstar/vecnorm(nstar);
            nstar = nmax^2*nstar/max(nstar); % 0 <= nstar <= nmax
            nstar(nstar<0) = 0; 
            nstar = sqrt(nstar);
            nstar = round(dec*nstar)/dec;
            alphastar = [atan2(Fdec(2),Fdec(1)) atan2(Fdec(4),Fdec(3)) atan2(Fdec(6),Fdec(5)) atan2(Fdec(8),Fdec(7))];
            alphastar = round(dec*alphastar)/dec;
            
            % calculate angular and linear acceleration with this alphastar and nstar
            [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alphastar, beta, teta,nstar, L, g, Mb, Mp, R);
            pdotdot = round(dec*pdotdot)/dec;
            
            M = [M Ib*wbdot];% Force produced by the MAV
            Meff = [Meff L*kf*nstar.'*nstar];
            
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

[Fmax, Ifmax] = max(vecnorm(F));
[Fmin, Ifmin] = min(vecnorm(F));
[Mmax, Immax] = max(vecnorm(M));
[Mmin, Immin] = min(vecnorm(M));

% Dfmin = D(:,Ifmin)
% Dfmax = D(:,Ifmax)
% Dmmin = D(:,Immin)
% Dmmax = D(:,Immax)

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

    % Generate a sphere consisting of 20by 20 faces of radius Fmin
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
% First Simple design representation in grey (beta = teta = 0)

% Generate a sphere
[x,y,z]=sphere;
% use surf function to plot
R= R/3;
r = 2*R/5;
propelerSphere1=surf(r*x+L,r*y,r*z); % centered at Op1
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

%The actual design with tilted arms
centerSphere=surf(R*x,R*y,R*z);% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
propelerSphere1=surf(r*x+Op1(1),r*y+Op1(2),r*z+Op1(3)); % centered at Op1
set(propelerSphere1,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere2=surf(r*x+Op2(1),r*y+Op2(2),r*z+Op2(3)); % centered at Op2
set(propelerSphere2,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere3=surf(r*x+Op3(1),r*y+Op3(2),r*z+Op3(3)); % centered at Op3
set(propelerSphere3,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere4=surf(r*x+Op4(1),r*y+Op4(2),r*z+Op4(3)); % centered at Op4
set(propelerSphere4,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
plot3([0 Op1(1)], [0 Op1(2)], [0 Op1(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op2(1)], [0 Op2(2)], [0 Op2(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op3(1)], [0 Op3(2)], [0 Op3(3)], 'c', 'LineWidth', 400*R)
plot3([0 Op4(1)], [0 Op4(2)], [0 Op4(3)], 'c', 'LineWidth', 400*R)
quiver3(0, 0, 0, 0, 0, -0.1, 'b')

x0=10;
y0=10;
width=1020;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height]);
daspect([1 1 1]);
title(['Design ' num2str(ii) ' (beta = [' num2str(rad2deg(beta(1))) ', ' num2str(rad2deg(beta(2))) ', ' num2str(rad2deg(beta(3))) ', ' num2str(rad2deg(beta(4))) '])'])
xlabel('X') 
ylabel('Y')
zlabel('Z')
camlight
end
end
