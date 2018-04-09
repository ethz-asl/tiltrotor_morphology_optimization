function [inv_q] = Quaternion_inverse(q)
%Quaternion_inverse Summary of this function goes here
%   Detailed explanation goes here
inv_q = [q(1) -q(2) -q(3) -q(4)]./(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
end

