function [Roty] = Roty(alpha)
%ROTY creates a rotation matrix, rotating alpha radians about the y-axis.

Roty = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];

end

