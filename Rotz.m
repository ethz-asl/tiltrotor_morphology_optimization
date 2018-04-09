function [RotZ] = RotZ(alpha)
%ROTZ Summary of this function goes here
%   Detailed explanation goes here
RotZ = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

end

