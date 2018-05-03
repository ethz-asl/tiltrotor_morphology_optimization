function [alpha, n, wRb, wb, p,t,P, Orientation] = Quadcopter_tilted_arms_simulation(start_time, end_time, dt, kf, km, wRb, wb, pdot, p, alpha, beta, theta,n, nmax, nmin, nhover, L, g, Mb, Mp, R, Rd, alphadotmax, dec, Op1, Op2, Op3, Op4)
%QUADCOPTER_TILTED_ARMS_SIMULATION Summary of this function goes here
%   Detailed explanation goes here
%% Loop to create a direction vector D with the desired number of direction:
% Simulation time [s]
times = start_time:dt:end_time;
nsquare = n.^2;
P = [];
Orientation = [];


%% Plot drone representation
figure(1);
[x,y,z]=sphere;
% use surf function to plot
R= R/3;
r = 2*R/5;
Op1 = wRb*(Op1+p);
Op2 = wRb*(Op2+p);
Op3 = wRb*(Op3+p);
Op4 = wRb*(Op4+p);
centerSphere=surf(R*x+p(1),R*y+ p(2),R*z+ p(3));% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
propelerSphere1=surf(r*x+Op1(1),r*y+Op1(2),r*z+Op1(3)); % centered at Op1
% if ~isequal(Op1, [L 0 0].')
%     % plot the original quadcopter design 
%     propelerSphere1=surf(r*x+L,r*y,r*z); hold on;% centered at Op1
%     set(propelerSphere1,'FaceColor',[0 0 0], ...
%        'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
%     plot3([0 L], [0 0], [0 0],'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
% end
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
% Plot drone arms
plot3([p(1) Op1(1)], [p(2) Op1(2)], [p(3) Op1(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op2(1)], [p(2) Op2(2)], [p(3) Op2(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op3(1)], [p(2) Op3(2)], [p(3) Op3(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op4(1)], [p(2) Op4(2)], [p(3) Op4(3)], 'c', 'LineWidth', 400*R)
% Plot drone propeller orientation
bRp1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
bRp2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
bRp3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
bRp4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
valpha1 = wRb*bRp1*[0; 0; 0.05];
valpha2 = wRb*bRp2*[0; 0; 0.05];
valpha3 = wRb*bRp3*[0; 0; 0.05];
valpha4 = wRb*bRp4*[0; 0; 0.05];
quiver3(Op1(1), Op1(2), Op1(3), valpha1(1), valpha1(2), valpha1(3), 'r')
quiver3(Op2(1), Op2(2), Op2(3), valpha2(1), valpha2(2), valpha2(3), 'r')
quiver3(Op3(1), Op3(2), Op3(3), valpha3(1), valpha3(2), valpha3(3), 'r')
quiver3(Op4(1), Op4(2), Op4(3), valpha4(1), valpha4(2), valpha4(3), 'r')

daspect([1 1 1]);
title('Illustration of the design')
xlabel('X') 
ylabel('Y')
zlabel('Z')
camlight
ii = 0;
% Step through the simulation, updating the state.
for t = times
    % Compute linear and angular accelerations.
    [m, Ib, pdotdot, wbdot, Op1, Op2, Op3, Op4] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta,n, L, g, Mb, Mp, R, true);
    wb = wb + dt * wbdot;
    wb = round(wb*dec)/dec;
    wb = mod(wb,2*pi);
    wRbdot = wRb*skew(wb);
    wRb = wRb + dt*wRbdot;
    wRb = round(wRb*dec)/dec;
    pdot = pdot + dt * pdotdot;
    pdot = round(pdot*dec)/dec;
    p = p + dt * pdot;
    p = round(p*dec)/dec;
    P = [P, p];
    [n,alphadot] = Quadcopter_tilted_arms_control(kf, km, L,m,Ib, wRb, wRbdot,-pdotdot, -pdot, -p, -wbdot, -wb, Rd, alpha,beta, theta, n, nmax, nmin, nhover, dt);
    n = round(dec*n)/dec;
    for i = 1:4
        if alphadot(i) >= alphadotmax 
            alphadot(i) = alphadotmax;
        elseif alphadot(i) <= -alphadotmax
             alphadot(i) = -alphadotmax;
        end
    end
    alphadot = round(alphadot*dec)/dec;
    alpha = alpha + dt * alphadot.';
    for i = 1:4
        if alpha(i) >= pi 
            alpha(i) = pi;
        elseif alpha(i) <= -pi
             alpha(i) = -pi;
        end
    end
    Orientation = wRb*[1 0 0].';% Orientation = [Orientation wRb*[1 0 0].'];
    if isequal(round(pdot*10^3)/10^3, [0;0;0]) && isequal(round(wb*10^3)/10^3, [0;0;0]) && t> 0
        break;
    end
    ii = ii+1;
    if ii == 10
        %% Plot drone representation
        figure(1);
        Op1 = wRb*(Op1+p);
        Op2 = wRb*(Op2+p);
        Op3 = wRb*(Op3+p);
        Op4 = wRb*(Op4+p);
        centerSphere=surf(r*x+p(1),r*y+ p(2),r*z+ p(3));% center of mass sphere
        set(centerSphere,'FaceColor',[0 0 0], ...
           'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none'); hold on;
        % Plot drone arms
        plot3([p(1) Op1(1)], [p(2) Op1(2)], [p(3) Op1(3)], 'k', 'LineWidth', 40*R)
        plot3([p(1) Op2(1)], [p(2) Op2(2)], [p(3) Op2(3)], 'k', 'LineWidth', 40*R)
        plot3([p(1) Op3(1)], [p(2) Op3(2)], [p(3) Op3(3)], 'k', 'LineWidth', 40*R)
        plot3([p(1) Op4(1)], [p(2) Op4(2)], [p(3) Op4(3)], 'k', 'LineWidth', 40*R)
        % Plot drone propeller orientation
        bRp1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
        bRp2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
        bRp3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
        bRp4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
        valpha1 = wRb*bRp1*[0; 0; 0.05];
        valpha2 = wRb*bRp2*[0; 0; 0.05];
        valpha3 = wRb*bRp3*[0; 0; 0.05];
        valpha4 = wRb*bRp4*[0; 0; 0.05];
        quiver3(Op1(1), Op1(2), Op1(3), valpha1(1), valpha1(2), valpha1(3), 'r')
        quiver3(Op2(1), Op2(2), Op2(3), valpha2(1), valpha2(2), valpha2(3), 'r')
        quiver3(Op3(1), Op3(2), Op3(3), valpha3(1), valpha3(2), valpha3(3), 'r')
        quiver3(Op4(1), Op4(2), Op4(3), valpha4(1), valpha4(2), valpha4(3), 'r')
        ii = 0;
    end
end
%% Plot drone representation
figure(1);
Op1 = wRb*(Op1+p);
Op2 = wRb*(Op2+p);
Op3 = wRb*(Op3+p);
Op4 = wRb*(Op4+p);
centerSphere=surf(R*x+p(1),R*y+ p(2),R*z+ p(3));% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
propelerSphere1=surf(r*x+Op1(1),r*y+Op1(2),r*z+Op1(3)); % centered at Op1
% if ~isequal(Op1, [L 0 0].')
%     % plot the original quadcopter design 
%     propelerSphere1=surf(r*x+L,r*y,r*z); hold on;% centered at Op1
%     set(propelerSphere1,'FaceColor',[0 0 0], ...
%        'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
%     plot3([0 L], [0 0], [0 0],'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
% end
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
% Plot drone arms
plot3([p(1) Op1(1)], [p(2) Op1(2)], [p(3) Op1(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op2(1)], [p(2) Op2(2)], [p(3) Op2(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op3(1)], [p(2) Op3(2)], [p(3) Op3(3)], 'c', 'LineWidth', 400*R)
plot3([p(1) Op4(1)], [p(2) Op4(2)], [p(3) Op4(3)], 'c', 'LineWidth', 400*R)
% Plot drone propeller orientation
bRp1 = rotz(rad2deg(theta(1)))*roty(rad2deg(beta(1)))*rotx(rad2deg(alpha(1)));
bRp2 = rotz(rad2deg(pi/2+theta(2)))*roty(rad2deg(beta(2)))*rotx(rad2deg(alpha(2)));
bRp3 = rotz(rad2deg(pi+theta(3)))*roty(rad2deg(beta(3)))*rotx(rad2deg(alpha(3)));
bRp4 = rotz(rad2deg(3*pi/2+theta(4)))*roty(rad2deg(beta(4)))*rotx(rad2deg(alpha(4)));
valpha1 = wRb*bRp1*[0; 0; 0.05];
valpha2 = wRb*bRp2*[0; 0; 0.05];
valpha3 = wRb*bRp3*[0; 0; 0.05];
valpha4 = wRb*bRp4*[0; 0; 0.05];
quiver3(Op1(1), Op1(2), Op1(3), valpha1(1), valpha1(2), valpha1(3), 'r')
quiver3(Op2(1), Op2(2), Op2(3), valpha2(1), valpha2(2), valpha2(3), 'r')
quiver3(Op3(1), Op3(2), Op3(3), valpha3(1), valpha3(2), valpha3(3), 'r')
quiver3(Op4(1), Op4(2), Op4(3), valpha4(1), valpha4(2), valpha4(3), 'r')
axis([-1 1 -1 1 0 3])
end

