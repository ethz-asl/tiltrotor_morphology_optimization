function [qtot] = Quaternion_product(q1,q2)
%Quaternion_product Hamilton product of two rotation quaternion
%The product of two rotation quaternions will be equivalent to the rotation
% q2 followed by the rotation q1
     
% vec = s1*v2 + s2*v1 + cross(v1,v2)
% scalar = s1*s2 - dot(v1,v2)
qtot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4); ...
        q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3);...
        q1(1)*q2(3)-q1(2)*q2(4)+q1(3)*q2(1)+q1(4)*q2(2); ...
        q1(1)*q2(4)+q1(2)*q2(3)-q1(3)*q2(2)+q1(4)*q2(1)];
end

