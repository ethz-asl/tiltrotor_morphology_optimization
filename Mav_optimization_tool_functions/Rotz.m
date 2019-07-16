function [RotZ] = RotZ(alpha)
%ROTZ creates a rotation matrix, rotating alpha radians about the z-axis.

RotZ = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

end

