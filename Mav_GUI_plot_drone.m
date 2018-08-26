function Mav_GUI_plot_drone(n, theta, beta, L, figure)
hold(figure,'off');
interval = 2*pi/n; % interval between arms in normal n-copter configuration
% pre-allocation:
beta = deg2rad(beta);
theta = deg2rad(theta);
bRp = zeros(3,3,n);
Op = zeros(3,n);
for i =1:n
    %% Find the propellers rotation matrix
    bRp(:,:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i));
    
    %% Find the propellers positions in the body frame
    Op(:,i) = Rotz((i-1)*interval)*Rotz(theta(i))*Roty(beta(i))*[L 0 0].';
end

wRb = eye(3);

%% Plot drone representation
% Generate a sphere
[x,y,z]=sphere;
% use surf function to plot
R = 0.06;
r = 2*R/5;
centerSphere=surf(figure, R*x,R*y,R*z);% center of mass sphere
set( centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold(figure,'on');
Op = wRb*Op;
Op0 = zeros(3,n);
interval = 2*pi/n;

for i = 1:n
    Op0(:,i) = wRb*Rotz((i-1)*interval)*[L 0 0].';
    propelerSphere=surf(figure, r*x+Op(1,i),r*y+Op(2,i),r*z+Op(3,i)); % centered at Op1
    set( propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
if ~isequal(Op(:,i), Op0(:,i))
    % plot the original quadcopter design 
    propelerSphere=surf(figure, r*x/2+Op0(1,i),r*y/2+Op0(2,i),r*z/2+Op0(3,i));% centered at Op0
    set( propelerSphere,'FaceColor',[.2 .2 .2], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')

    
    plot3(figure, [0 Op0(1,i)], [0 Op0(2,i)], [0 Op0(3,i)],'Color',[0.5 0.5 0.5], 'LineWidth', 10*R)
end
set(propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')

plot3(figure, [0 Op(1,i)], [0 Op(2,i)], [0 Op(3,i)], 'c', 'LineWidth', 100*R)

%% Plot thruster direction
TD = wRb*bRp(:,:,i)*[0; 0; L/5];
text(figure, Op(1,i)+TD(1), Op(2,i)+TD(2),Op(3,i) + TD(3), ['p_{' num2str(i) '}']) 


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
    angletheta = fill3(figure, theta0(1,:),theta0(2,:),theta0(3,:),'r');
    set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op1+Op00)/2;
    text(figure, postxt(1), postxt(2), postxt(3), ['\theta_{' num2str(i) '}'])
end

%% plot angle beta
if beta(i) ~= 0
    Op00 = Rotz((i-1)*interval)*Rotz(theta(i))*[L 0 0].';
    Op00 = cos(beta(i))*Op00;
    beta0 = [[0; 0; 0], Op00, Op(:,i) ];
    anglebeta = fill3(figure, beta0(1,:),beta0(2,:),beta0(3,:),'b');
    set(anglebeta,'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op(:,i)+Op00)/2;
    text(figure, postxt(1), postxt(2), postxt(3), ['\beta_{' num2str(i) '}']) 
end

end

%% plot drone axis
quiver3(figure, 0, 0, 0, wRb(1,1)*0.25, wRb(2,1)*0.25,wRb(3,1)*0.25, 'k')
text(figure, wRb(1,1)*0.25, wRb(2,1)*0.25, wRb(3,1)*0.25, 'x') 
quiver3(figure, 0, 0, 0, wRb(1,2)*0.25, wRb(2,2)*0.25, wRb(3,2)*0.25, 'k')
text(figure, wRb(1,2)*0.25, wRb(2,2)*0.25, wRb(3,2)*0.25, 'y')
quiver3(figure, 0, 0, 0, wRb(1,3)*0.25, wRb(2,3)*0.25, wRb(3,3)*0.25, 'k')
text(figure, wRb(1,3)*0.25, wRb(2,3)*0.25, wRb(3,3)*0.25, 'z')


daspect(figure, [1 1 1]);
title(figure, '')
text(figure, -L , 0, L+0.2,'Initial solution representation:')
end