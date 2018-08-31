function Mav_plot(n, wRb, fig_number, design_number, theta, beta,  D_unit, F, Feff, M,Meff, Heff, L, Op, bRp, worthF, worthM, worthH, number_of_directions, plot_volume, TRI, F_surf, F_vol, M_surf, M_vol)
%MAV_PLOT Summary of this function goes here
%   Detailed explanation goes here
figure(fig_number); 
%% Plot Force space
subplot(2,2,1);
if  plot_volume
    % Surface Reconstruction from scattered points cloud
    trisurf(TRI,F(1,:),F(2,:),F(3,:),'facecolor','k', 'FaceAlpha', 0.3, 'edgecolor','none'); hold on;
    scatter_size = 50;
else
    scatter_size = 100;
end
colormap(flipud(jet));
scatter3(F(1,:), F(2,:), F(3,:),  scatter_size ,Feff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of the Thrust (%)';
caxis([min(Feff) max(Feff)]);
daspect([1 1 1]);
title('Reachable force space', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% Plot torque space
subplot(2,2,2);
if  plot_volume
    % Surface Reconstruction from scattered points cloud
    trisurf(TRI,M(1,:),M(2,:),M(3,:),'facecolor','k', 'FaceAlpha', 0.3, 'edgecolor','none'); hold on;
    scatter_size = 50;
else
    scatter_size = 100;
end
colormap(flipud(jet));
scatter3(M(1,:), M(2,:), M(3,:),  scatter_size ,Meff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of the Torque (%)';
caxis([min(Meff) max(Meff)])
daspect([1 1 1]);
title('Reachable torque space', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% Plot Hover efficincy
subplot(2,2,3);
if  plot_volume
    % Surface Reconstruction from scattered points cloud
    trisurf(TRI,D_unit(1,:),D_unit(2,:),D_unit(3,:),'facecolor','k', 'FaceAlpha', 0.3, 'edgecolor','none'); hold on;
    scatter_size = 50;
else
    scatter_size = 100;
end
colors = flipud(jet);
colormap(colors);
scatter3(D_unit(1,:), D_unit(2,:), D_unit(3,:),  scatter_size ,Heff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of Hover (%)';
caxis([min(Heff) max(Heff)])
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
R = L*0.1/0.5;
r = 2*R/5;
centerSphere=surf(R*x,R*y,R*z);% center of mass sphere
set(centerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none'); hold on;
Op = wRb*Op;
Op0 = zeros(3,n);
interval = 2*pi/n;

for i = 1:n
    Op0(:,i) = wRb*Rotz((i-1)*interval)*[L 0 0].';
    propelerSphere=surf(r*x+Op(1,i),r*y+Op(2,i),r*z+Op(3,i)); % centered at Op1
    set(propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
if ~isequal(Op(:,i), Op0(:,i))
    % plot the original quadcopter design 
    propelerSphere=surf(r*x/2+Op0(1,i),r*y/2+Op0(2,i),r*z/2+Op0(3,i)); hold on;% centered at Op0
    set(propelerSphere,'FaceColor',[.2 .2 .2], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')

    
    plot3([0 Op0(1,i)], [0 Op0(2,i)], [0 Op0(3,i)],'Color',[0.5 0.5 0.5], 'LineWidth', 10*R)
end
set(propelerSphere,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')

plot3([0 Op(1,i)], [0 Op(2,i)], [0 Op(3,i)], 'c', 'LineWidth', 100*R)

%% Plot thruster direction
TD = wRb*bRp(:,:,i)*[0; 0; L/5];
text( Op(1,i)+TD(1), Op(2,i)+TD(2),Op(3,i) + TD(3), ['P_{' num2str(i) '}']) 


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
    text(postxt(1), postxt(2), postxt(3), ['\theta_{' num2str(i) '}'])
end

%% plot angle beta
if beta(i) ~= 0
    Op00 = Rotz((i-1)*interval)*Rotz(theta(i))*[L 0 0].';
    Op00 = cos(beta(i))*Op00;
    beta0 = [[0; 0; 0], Op00, Op(:,i) ];
    fill3(beta0(1,:),beta0(2,:),beta0(3,:),'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    postxt = (Op(:,i)+Op00)/2;
    text(postxt(1), postxt(2), postxt(3), ['\beta_{' num2str(i) '}']) 
end

end

%% plot drone axis
quiver3(0, 0, 0, wRb(1,1)*0.1, wRb(2,1)*0.1,wRb(3,1)*0.1, 'k')
text(wRb(1,1)*0.1, wRb(2,1)*0.1, wRb(3,1)*0.1, 'x') 
quiver3(0, 0, 0, wRb(1,2)*0.1, wRb(2,2)*0.1, wRb(3,2)*0.1, 'k')
text(wRb(1,2)*0.1, wRb(2,2)*0.1, wRb(3,2)*0.1, 'y')
quiver3(0, 0, 0, wRb(1,3)*0.1, wRb(2,3)*0.1, wRb(3,3)*0.1, 'k')
text(wRb(1,3)*0.1, wRb(2,3)*0.1, wRb(3,3)*0.1, 'z')

%% Plot gravity direction

daspect([1 1 1]);
title('Illustration of the design')
xlabel('X')
ylabel('Y')
zlabel('Z')
camlight

%% General plot options and annotations
x0=10;
y0=10;
width=1800;
height= 900;
set(gcf,'units','points','position',[x0,y0,width,height]);

Fnorm = vecnorm(F);
Fmean = mean(Fnorm);
Fmad = mad(Fnorm);
Mnorm = vecnorm(M);
Mmean = mean(Mnorm);
Mmad = mad(Mnorm);
Hmean = mean(Heff);
Hmad = mad(Heff);
str = (['Design ' num2str(design_number) ': \beta = ' mat2str(round(rad2deg(beta)*10^2)/10^2) ...
    ', \theta = ' mat2str(round(rad2deg(theta)*10^2)/10^2) ', L = ' num2str(round(10^2*L)/10^2)...
    '. The reacheable force space: surface = ' num2str(round(10^2*F_surf)/10^2) ' [N^{2}]'  ...
    ', volume = ' num2str(round(10^2*F_vol)/10^2) ' [N^{3}], mean = '  num2str(round(10^2*Fmean)/10^2) ...
    ' [N] and mean absolute deviation = ' num2str(round(10^2*Fmad)/10^2) '[N]' ...
    ', with F_{min} = ' num2str(min(vecnorm(F))) ', F_{max} = ' num2str(max(vecnorm(F))) ...
    '. The reacheable torque space: surface = ' num2str(round(10^2*M_surf)/10^2) ' [Nm^{2}]'  ...
    ', volume = ' num2str(round(10^2*M_vol)/10^2) ' [Nm^{3}], mean = '  num2str(round(10^2*Mmean)/10^2) ...
    ' [Nm] and mean absolute deviation = ' num2str(round(10^2*Mmad)/10^2) '[Nm]' ...
    ', with M_{min} = ' num2str(min(vecnorm(M))) ', M_{max} = ' num2str(max(vecnorm(M))) ...
    '. The hover efficiency, mean = '  num2str(round(10^2*Hmean)/10^2) '[%]' ...
    ' and mean absolute deviation = ' num2str(round(10^2*Hmad)/10^2) '[%]' ...
    ', with H_{min} = ' num2str(min(Heff)) ', H_{max} = ' num2str(max(Heff))]);
dim = [ .1 0.955 .9 .045];
annotation('textbox',dim,'String',str,'FitBoxToText','off');
if worthF ~= 0 || worthM ~= 0 || worthH ~= 0
    str = (['The optimization improved the maximum force in ' num2str(worthF) ' directions'  ...
        ', the maximum torque in ' num2str(worthM) ' directions'  ...
        ', the efficiency of hover in ' num2str(worthH) ' directions' ...
        ', on a total of ' num2str(number_of_directions) ' directions']);
    dim = [ .3 .5 .9 .025];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
end
end

