clear all;
close all;

[file_path] = fileparts(mfilename('fullpath'));
addpath(file_path);
addpath([file_path '/Mav_optimization_tool_functions/']);


figure(1);
n = 6;
L = 0.5;

theta = zeros(1,n);
beta = zeros(1,n);

beta = ones(1,n) * deg2rad(30);
theta = ones(1,n) * deg2rad(20);

% beta = deg2rad([0 0 36.89 -58.4 41.45]);
% theta = deg2rad([0 30.3 -16.13 -5.46 -33.37]);

% beta = -deg2rad([32.42, 35.49, 35.44, 35.49]);

% beta = acos(sqrt(2/3))*[-1, -1, -1, -1, 1, 1, 1, 1];
% theta = [0, pi/4, pi/2, 3*pi/4, -pi, -3*pi/4, -pi/2, -pi/4];
alpha = zeros(1,n);
w = zeros(1,n);
wRb = eye(3);
[g, dec, kf, km] = Mav_parameters();
[m, Ib, pdotdot, wbdot, Op, bRp] = Mav_dynamic(n, kf, km, wRb, alpha, beta, theta,w, L, g, dec, false);
%% Plot drone representation
% Generate a sphere
[x,y,z]=sphere;
% use surf function to plot
R = L*0.1/0.6;
r = 2*R/5;
centerSphere=surf(R*x,R*y,R*z);% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',01.0,'FaceLighting','gouraud','EdgeColor','none'); hold on;
Op = wRb*Op;
Op0 = zeros(3,n);
interval = 2*pi/n;

for i = 1:n
    Op0(:,i) = wRb*Rotz((i-1)*interval)*[L 0 0].';
    propelerSphere=surf(r*x+Op(1,i),r*y+Op(2,i),r*z+Op(3,i)); % centered at Op1
    set(propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')

set(propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')

plot3([0 Op(1,i)], [0 Op(2,i)], [0 Op(3,i)], 'c', 'LineWidth', 100*R)

%% Plot propeller number
if mod(n,2) == 0
   if mod(i,2) == 0
    TD = wRb*bRp(:,:,i)*[0; 0; L/5];
   else
   	TD = bRp(:,:,i)*[0; 0; -L/5];
   end
else
    TD = bRp(:,:,i)*[0; 0; -L/5];
end

% if i<5
%     TD = wRb*bRp(:,:,i)*[0; 0; L/5];
% else
%     TD = bRp(:,:,i)*[0; 0; -L/5];
% end
% T = text( Op(1,i)+TD(1), Op(2,i)+TD(2),Op(3,i) + TD(3), ['R_{' num2str(i) '}']); 
% set(T, ...
%     'FontName'   , 'Modern No. 20' , ...
%     'FontWeight' , 'normal'          , ...
%     'FontSize'   , 20        );

%% plot angle theta
if theta(i) ~= 0
    Op00 = Op0(:,i); % Op for a normal design of n-copter
    Op00norm = norm(Op00);
    Op00(3) = 0;
    Op00 = cos(beta(i))*Op00*Op00norm/norm(Op00);
    
    Op1 = Op(:,i);
    Op1norm = norm(Op1);
    Op1(3) = 0;
    Op1 = cos(beta(i))*Op1*Op1norm/norm(Op1);
    
    theta0 = [[0; 0; 0], Op00, Op1 ];
    fill3(theta0(1,:),theta0(2,:),theta0(3,:),'r', 'FaceAlpha', 0.2, 'EdgeColor','none'); hold on;
    postxt = (Op1+Op00)/2;
%     T = text(postxt(1), postxt(2), postxt(3), ['\theta_{' num2str(i) '}']);
%     set(T, ...
%     'FontName'   , 'Modern No. 20' , ...
%     'FontWeight' , 'normal'          , ...
%     'FontSize'   , 30        );
end

%% plot angle beta
if beta(i) ~= 0
    Op00 = Rotz((i-1)*interval)*Rotz(theta(i))*[L 0 0].';
    Op00 = cos(beta(i))*Op00;
    beta0 = [[0; 0; 0], Op00, Op(:,i) ];
    fill3(beta0(1,:),beta0(2,:),beta0(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op(:,i)+Op00)/2;
%     T = text(postxt(1), postxt(2), postxt(3), ['\beta_{' num2str(i) '}']);
%     set(T, ...
%     'FontName'   , 'Modern No. 20' , ...
%     'FontWeight' , 'normal'          , ...
%     'FontSize'   , 30        );
end

end

%% plot drone axis
quiver3(-(L), 0, 0, 2*(L+0.075), 0, 0, 'r', 'LineWidth', 30*R)
T1 = text((L+0.1), 0, 0, 'X');
set(T1, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'normal'          , ...
    'FontSize'   , 30        );
quiver3(0, -(L), 0, 0, 2*(L+0.075), 0,  'g', 'LineWidth', 30*R)
T2 = text(0, (L+0.1), 0, 'Y');
set(T2, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'normal'          , ...
    'FontSize'   , 30       );
quiver3(0, 0, 0, 0, 0, (L), 'b', 'LineWidth', 30*R)
T3 = text(0, 0, (L), 'Z');
set(T3, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'normal'          , ...
    'FontSize'   , 30        );
%% Plot platonic solid
% Op = [Op, Op(:,1), Op(:,3), Op(:,2), Op(:,4)]; % QUADCOPTER
% Op = [Op, Op(:,1), Op(:,3), Op(:,5), Op(:,1), Op(:,6), Op(:,4), Op(:,2), Op(:,6)]; % HEXACOPTER
% Op = [Op(:,1), Op(:,2), Op(:,6), Op(:,5) Op(:,1), Op(:,4), Op(:,8), Op(:,5), Op(:,8), Op(:,7), Op(:,3), Op(:,4), Op(:,3), Op(:,2), Op(:,6), Op(:,7)]; % OCTACOPTER
% patch(Op(1,:), Op(2,:), Op(3,:),'k', 'FaceAlpha', 0.0, 'EdgeColor','k');

% patch(0.4*[1 -1 -1 1], 0.4*[1 1 -1 -1], 0.4*[0 0 0 0],'k', 'FaceAlpha', 0.4, 'EdgeColor','none')  
daspect([1 1 1]);
camlight

%% General plot options and annotations
x0=10;
y0=10;
width=1800;
height= 900;
set(gcf,'units','points','position',[x0,y0,width,height]);
hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none';  
axis off