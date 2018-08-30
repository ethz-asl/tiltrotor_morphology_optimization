close all
clear all
figure(1)
z = 1/sqrt(2);
S.Vertices = ([1,0,-z;-1,0,-z;0,1,z;0,-1,z]')';
S.Faces = [1,3,4;2,3,4;1,2,3;1,2,4];
S.FaceVertexCData = [ 1 ];
p = patch(S); hold on
set( p,'FaceColor',[1 1 1], ...
    'FaceAlpha',0.0,'FaceLighting','gouraud','EdgeColor','k')

[X,Y,Z]=sphere;
% use surf function to plot
R = 0.05;
centerSphere=surf(R*X,R*Y,R*Z);% center of mass sphere
set( centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');




theta0 = [[0; 0; 0], [1; 0;-z], [1; 0; 0]];
angletheta = fill3(theta0(1,:),theta0(2,:),theta0(3,:),'b');
set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
postxt = (([1; 0;-z]+ [1; 0; 0])/2);
Text = text(postxt(1), postxt(2), postxt(3), ['\beta'])

set(Text, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'bold'      );

set( Text                    , ...
    'FontSize'   , 30         , ...
    'FontWeight' , 'normal'   );

theta0 = [[0; 0; 0], [0; 1;z], [0; 1; 0]];
angletheta = fill3(theta0(1,:),theta0(2,:),theta0(3,:),'b');
set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
postxt = ([0,1,z] + [0, 1, 0])/2;
Text = text(postxt(1), postxt(2), postxt(3), ['\beta'])

set(Text, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'bold'      );

set( Text                    , ...
    'FontSize'   , 30         , ...
    'FontWeight' , 'normal'   );

% set(figure(1), 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')

%%cube
figure(2)
xc=0; yc=0; zc=0;    % coordinated of the center
L=1;                 % cube size (length of an edge)


X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

C='blue';                  % unicolor

X = L*(X-0.5) + xc;
Y = L*(Y-0.5) + yc;
Z = L*(Z-0.5) + zc; 

p = fill3(X,Y,Z,C); hold on   % draw cube

set( p,'FaceColor',[1 1 1], ...
    'FaceAlpha',0.0,'FaceLighting','gouraud','EdgeColor','k')


theta0 = [[0; 0; 0], [0.5;0.5;0], [0.5;0.5;-0.5]];
angletheta = fill3(theta0(1,:),theta0(2,:),theta0(3,:),'b');
set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
postxt = (([0.5;0.5;0]+ [0.5;0.5;-0.5])/2);
Text = text(postxt(1), postxt(2), postxt(3), ['\beta'])

set(Text, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'bold'      );

set( Text                    , ...
    'FontSize'   , 30         , ...
    'FontWeight' , 'normal'   );

theta0 = [[0; 0; 0], [-0.5;0.5;0], [-0.5;0.5;0.5]];
angletheta = fill3(theta0(1,:),theta0(2,:),theta0(3,:),'b');
set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
postxt = ([-0.5;0.5;0]+ [-0.5;0.5;0.5])/2;
Text = text(postxt(1), postxt(2), postxt(3), ['\beta'])

set(Text, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'bold'      );

set( Text                    , ...
    'FontSize'   , 30         , ...
    'FontWeight' , 'normal'   );
[X,Y,Z]=sphere;

% use surf function to plot
R = 0.05;
centerSphere=surf(R*X,R*Y,R*Z);% center of mass sphere
set( centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');

figure(3)

[V,F]=platonic_solid(3,1);
patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.0,'EdgeColor','k','LineWidth',1);
hold on
[X,Y,Z]=sphere;

% use surf function to plot
R = 0.05;
centerSphere=surf(R*X,R*Y,R*Z);% center of mass sphere
set( centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');


theta0 = [[0; 0; 0], [-0.707;0.707;0], [-0.3535;0.3535;0.5]];
angletheta = fill3(theta0(1,:),theta0(2,:),theta0(3,:),'b');
set(angletheta,'FaceAlpha', 0.2, 'EdgeColor','none');
postxt = ([-0.5;0.5;0]+ [-0.5;0.5;0.5])/2;
Text = text(postxt(1), postxt(2), postxt(3), ['\beta'])

set(Text, ...
    'FontName'   , 'Modern No. 20' , ...
    'FontWeight' , 'bold'      );

set( Text                    , ...
    'FontSize'   , 30         , ...
    'FontWeight' , 'normal'   );