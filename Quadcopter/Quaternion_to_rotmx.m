function [Rot] = Quaternion_to_rotmx(q)
%Quaternion_to_rotmx Summary of this function goes here
%   Detailed explanation goes here
% Reshape the quaternions in the depth dimension

s = q(1);
x = q(2);
y = q(3);
z = q(4);

% Explicitly define concatenation dimension for codegen
Rot = [1-2*(y^2 + z^2), 2*(x*y - s*z), 2*(x*z + s*y); ...
       2*(x*y + s*z), 1-2*(x^2 + z^2), 2*(y*z - s*x);...
       2*(x*z - s*y), 2*(y*z + s*x), 1-2*(x^2 + y^2)];

end

