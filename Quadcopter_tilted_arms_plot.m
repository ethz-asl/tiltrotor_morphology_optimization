function Quadcopter_tilted_arms_plot(fig_number, design_number, theta, beta, D, F, Feff, M,Meff, Heff, L, R, Op1, Op2, Op3, Op4, step, worthF, worthM, worthH, errF, errM, errH, ii)
%QUADCOPTER_TILTED_ARMS_PLOT Summary of this function goes here
%   Detailed explanation goes here
figure(fig_number); 
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
D_unit = D./vecnorm(D);
scatter3(D_unit(1,:), D_unit(2,:), D_unit(3,:),  100 ,Heff, 'filled'); hold on;
c = colorbar('eastoutside');
c.Label.String = 'Efficiency of Hover (%)';
caxis([min(Heff) max(Heff)])
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
    plot3([0 L], [0 0], [0 0],'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
end
set(propelerSphere1,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere2=surf(r*x+Op2(1),r*y+Op2(2),r*z+Op2(3)); % centered at Op2
if ~isequal(Op1, [0 L 0].')
    % plot the original quadcopter design 
    propelerSphere2=surf(r*x+0,r*y+L,r*z+0); % centered at Op2
    set(propelerSphere2,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 0], [0 L], [0 0], 'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
end
set(propelerSphere2,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere3=surf(r*x+Op3(1),r*y+Op3(2),r*z+Op3(3)); % centered at Op3
if ~isequal(Op1, [-L 0 0].')
    % plot the original quadcopter design 
    propelerSphere3=surf(r*x-L,r*y,r*z); % centered at Op3
    set(propelerSphere3,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 -L], [0 0], [0 0], 'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
end
set(propelerSphere3,'FaceColor',[0 0 0], ...
   'FaceAlpha',0.6,'FaceLighting','gouraud','EdgeColor','none')
propelerSphere4=surf(r*x+Op4(1),r*y+Op4(2),r*z+Op4(3)); % centered at Op4
if ~isequal(Op1, [0 -L 0].')
    % plot the original quadcopter design 
    propelerSphere4=surf(r*x,r*y-L,r*z); % centered at Op4
    set(propelerSphere4,'FaceColor',[0 0 0], ...
       'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
    plot3([0 0], [0 -L], [0 0], 'Color',[0.5 0.5 0.5], 'LineWidth', 40*R)
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
if theta(4) ~= 0
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
if beta(4) ~= 0
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

str = (['Design ' num2str(design_number) ': \beta = [' num2str(rad2deg(beta(1))) ', ' ...
    num2str(rad2deg(beta(2))) ', ' num2str(rad2deg(beta(3))) ...
    ', ' num2str(rad2deg(beta(4))) '], \theta = [' ...
    num2str(round(rad2deg(theta(1)))) ', ' num2str(round(rad2deg(theta(2)))) ', ' ...
    num2str(round(rad2deg(theta(3)))) ', ' num2str(round(rad2deg(theta(4)))) ...
    '], ' 'F_{min} = ' num2str(min(vecnorm(F))) ', F_{max} = ' num2str(max(vecnorm(F))) ...
    ', M_{min} = ' num2str(min(vecnorm(M))) ', M_{max} = ' num2str(max(vecnorm(M))) ...
    ', H_{min} = ' num2str(min(Heff)) ', H_{max} = ' num2str(max(Heff))]);
dim = [ .2 .7 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
str = (['Errors Thrust Opt: ' num2str(errF) ', Errors Torque Opt: ' ...
            num2str(errM) ', Errors Hover Opt: ' num2str(errH) ...
            ', Worth Opt Thrust : ' num2str(worthF) ', Worth Opt Torque: ' ...
            num2str(worthM) ', Worth Opt Hover: ' num2str(worthH) ', ii: ' num2str(ii)]);  
dim = [ .2 .25 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
end