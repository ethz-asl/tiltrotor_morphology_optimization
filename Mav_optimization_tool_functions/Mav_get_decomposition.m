function [w0,alpha0] = Mav_get_decomposition(n, dec, kf, Fdec)
%MAV_GET_DECOMPOSITION finds the rotors speeds and orientations from the
%decomposed force vector Fdec

% Inverse substitution :
    %                       wi² = sqrt(Fdec(2*i-1)² + Fdec(2*i)²)/kf  (Rotor speed [tour/s])
    %                       alphai = atan2(Fdec(2*i), Fdec(2*i-1))  (tilting angles [rad])
    
w0 = zeros(n,1);
alpha0 = zeros(n,1);
for i = 1:n
    w0(i) = 1/kf*sqrt(Fdec(2*i-1)^2 + Fdec(2*i)^2);
    alpha0(i) = atan2(Fdec(2*i),Fdec(2*i-1));
end
w0 = sqrt(w0);
w0 = round(dec*w0)/dec; % rotors orientations
alpha0 = round(dec*alpha0)/dec; % rotors speeds
end

