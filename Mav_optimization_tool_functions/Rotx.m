function [Rotx] = Rotx(alpha)
%ROTX creates a rotation matrix, rotating alpha radians about the x-axis.

Rotx = [1 0 0;0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];

end