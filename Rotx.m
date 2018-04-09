function [Rotx] = Rotx(alpha)
%ROTX Summary of this function goes here
%   Detailed explanation goes here

Rotx = [1 0 0;0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];

end