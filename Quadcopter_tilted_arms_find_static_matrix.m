%% Parameters
Mb = 1; % mass body [kg]
R = 0.05; % Radius of the body (Body assumed to be a sphere)
Mp = 0.1; % propeller mass [kg]
m = Mb+4*Mp;% drone mass [kg]
g = 9.81;
 kf = 3.86e-4; % Propeller thrust coefficient % [kg.m]
% km = 2e-5;% Propeller drag coefficient
% L = 0.15; % Arm length [m]
nmax =150; % [roun/s]
nhover = sqrt((m*g/4)/kf); % [roun/s]
Ndecimals = 4;
dec = 10.^Ndecimals;
syms kf km L;
%% Simulation
syms theta1 theta2 theta3 theta4;
theta = [theta1 theta2 theta3 theta4];
syms beta1 beta2 beta3 beta4;
beta = [beta1 beta2 beta3 beta4];
syms alpha1 alpha2 alpha3 alpha4;
alpha = [alpha1 alpha2 alpha3 alpha4];
syms n1 n2 n3 n4;
n = [n1 n2 n3 n4];

beta0 = [pi/6 -pi/6 -pi/6 pi/6];
theta0 = [0 0 0 0];
alpha0 = [0 0 0 0];
use_quaternions = false;

if ~use_quaternions
    %% Find the propellers frame rotation matrix
    bRp1 = Rotz(theta(1))*Roty(beta(1))*Rotx(alpha(1));
    bRp2 = Rotz(pi/2)*Rotz(theta(2))*Roty(beta(2))*Rotx(alpha(2));
    bRp3 = Rotz(pi)*Rotz(theta(3))*Roty(beta(3))*Rotx(alpha(3));
    bRp4 = Rotz(3*pi/2)*Rotz(theta(4))*Roty(beta(4))*Rotx(alpha(4));
    
    %% Find the propellers positions in the body frame
    Op1 = Rotz(theta(1))*Roty(beta(1))*[L 0 0].';
    Op2 = Rotz(pi/2)*Rotz(theta(2))*Roty(beta(2))*[L 0 0].';
    Op3 = Rotz(pi)*Rotz(theta(3))*Roty(beta(3))*[L 0 0].';
    Op4 = Rotz(3*pi/2)*Rotz(theta(4))*Roty(beta(4))*[L 0 0].';
    
    %% Find forces applied by all propellers thrusts on the body (in body frame)
    Tp1 = [0 0 kf*n1^2].';
    Tp2 = [0 0 kf*n2^2].';
    Tp3 = [0 0 kf*n3^2].';
    Tp4 = [0 0 kf*n4^2].';

    F  = bRp1*Tp1 + bRp2*Tp2 + bRp3*Tp3 + bRp4*Tp4


    %% Find torques applied by all propellers on the body (in body frame)
    Tauext1 = [0 0 -km*n(1)^2].'; % Thrust vector propeller 1
    Tauext2 = [0 0 km*n(2)^2].'; % Thrust vector propeller 2
    Tauext3 = [0 0 -km*n(3)^2].'; % Thrust vector propeller 3
    Tauext4 = [0 0 km*n(4)^2].'; % Thrust vector propeller 4
    % neglect rotor counter torque:
    % Tauext1 = [0 0 0].';
    % Tauext2 = [0 0 0].';
    % Tauext3 = [0 0 0].';
    % Tauext4 = [0 0 0].';
    Taub = cross(Op1,bRp1*Tp1) + cross(Op2,bRp2*Tp2) + cross(Op3,bRp3*Tp3) + cross(Op4,bRp4*Tp4);
    M = bRp1*Tauext1 + bRp2*Tauext2 + bRp3*Tauext3 + bRp4*Tauext4 + Taub;
        
    %% Find the force static matix 
    % matrix that link the force applied on the body to the vector:
    % kf*[cos(alpha1)n1^2; sin(alpha1)n1^2; cos(alpha2)n2^2; sin(alpha2)n2^2;
    %     cos(alpha3)n3^2; sin(alpha3)n3^2; cos(alpha4)n4^2; sin(alpha4)n4^2]
    
    [row columns] = size(n);
    A = [];
    Astatic = [];
    for ii = 1:3
        for jj = 1:columns
            Aij = F(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , ...
                {1,1,0}); Aij = subs(Aij,n , [0 0 0 0]);
            Ai2j = F(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)}, ...
                {1,1,pi/2}); Ai2j = subs(Ai2j,n , [0 0 0 0]);
            A = [A, Aij, Ai2j];
        end
         Astatic = [Astatic; A];
         A = [];
    end
    
    %% Find the torque static matix 
    % matrix that link the torque applied on the body to the vector:
    % kf*[cos(alpha1)n1^2; sin(alpha1)n1^2; cos(alpha2)n2^2; sin(alpha2)n2^2;
    %     cos(alpha3)n3^2; sin(alpha3)n3^2; cos(alpha4)n4^2; sin(alpha4)n4^2]
    
    A = [];
    A_M_stat = [];
    for ii = 1:3
        for jj = 1:columns
            Aij = M(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , {1,1,0});
            Aij = subs(Aij,n , [0 0 0 0]); Aij = subs(Aij,km , km/kf);
            Ai2j = M(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)} , {1,1,pi/2});
            Ai2j = subs(Ai2j,n , [0 0 0 0]); Ai2j = subs(Ai2j,km , km/kf);
            A = [A, Aij, Ai2j];
        end
         A_M_stat = [A_M_stat; A];
         A = [];
    end
    %% Test to verify solution
    
%     bRp1 = subs(bRp1, beta, beta0);
%     bRp1 = subs(bRp1, theta, theta0)
%     bRp2 = subs(bRp2, beta, beta0);
%     bRp2 = subs(bRp2, theta, theta0)
%     bRp3 = subs(bRp3, beta, beta0);
%     bRp3 = subs(bRp3, theta, theta0)
%     bRp4 = subs(bRp4, beta, beta0);
%     bRp4 = subs(bRp4, theta, theta0)
%     F = subs(F, beta, beta0);
%     F = subs(F, theta, theta0);
%     F = subs(F, alpha, alpha0)
%     Astatic1 = subs(Astatic, beta, beta0);
%     Astatic1 = subs(Astatic1, theta, theta0)
%     Astatic1 = double(round(Astatic1*10^5)/10^5)
    
%     M = subs(M, beta, beta0);
%     M = subs(M, theta, theta0);
%     M = subs(M, alpha, alpha0)
%     A_M_stat = subs(A_M_stat, beta, beta0);
%     A_M_stat = subs(A_M_stat, theta, theta0)

%     Op1 = subs(Op1, beta, beta0);
%     Op1 = subs(Op1, theta, theta0);
%     Op1 = subs(Op1, L, 1)
%     Op2 = subs(Op2, beta, beta0);
%     Op2 = subs(Op2, theta, theta0);
%     Op2 = subs(Op2, L, 1)
%     Op3 = subs(Op3, beta, beta0);
%     Op3 = subs(Op3, theta, theta0);
%     Op3 = subs(Op3, L, 1)
%     Op4 = subs(Op4, beta, beta0);
%     Op4 = subs(Op4, theta, theta0);
%     Op4 = subs(Op4, L, 1)
%     
%     Op1 = round(Op1*10^5)/10^5
%     Op2 = round(Op2*10^5)/10^5
%     Op3 = round(Op3*10^5)/10^5
%     Op4 = round(Op4*10^5)/10^5
    
else
    %% Find the propellers frame rotation quaternion
    qalpha1 = [cos(alpha(1)/2); sin(alpha(1)/2); 0; 0];
    qalpha2 = [cos(alpha(2)/2); sin(alpha(2)/2); 0; 0];
    qalpha3 = [cos(alpha(3)/2); sin(alpha(3)/2); 0; 0]; 
    qalpha4 = [cos(alpha(4)/2); sin(alpha(4)/2); 0; 0];

    qbeta1 = [cos(beta(1)/2); 0; sin(beta(1)/2); 0];
    qbeta2 = [cos(beta(2)/2); 0; sin(beta(2)/2); 0];
    qbeta3 = [cos(beta(3)/2); 0; sin(beta(3)/2); 0];
    qbeta4 = [cos(beta(4)/2); 0; sin(beta(4)/2); 0];

    q1 = Quaternion_product(qbeta1,qalpha1);
    q2 = Quaternion_product(qbeta2,qalpha2);
    q3 = Quaternion_product(qbeta3,qalpha3);
    q4 = Quaternion_product(qbeta4,qalpha4);

    qtheta1 = [cos(theta(1)/2); 0; 0; sin(theta(1)/2)];
    qtheta2 = [cos((pi/2+theta(2))/2); 0; 0; sin((pi/2+theta(2))/2)];
    qtheta3 = [cos((pi+theta(3))/2); 0; 0; sin((pi+theta(3))/2)];
    qtheta4 = [cos((3*pi/2+theta(4))/2); 0; 0; sin((3*pi/2+theta(4))/2)];

    q1 = Quaternion_product(qtheta1,q1);
    q2 = Quaternion_product(qtheta2,q2);
    q3 = Quaternion_product(qtheta3,q3);
    q4 = Quaternion_product(qtheta4,q4);

    inv_q1 = Quaternion_inverse(q1);
    inv_q2 = Quaternion_inverse(q2);
    inv_q3 = Quaternion_inverse(q3);
    inv_q4 = Quaternion_inverse(q4);

    %% Find the propellers positions in the body frame
    q01 = Quaternion_product(qtheta1, qbeta1);
    q02 = Quaternion_product(qtheta2, qbeta2);
    q03 = Quaternion_product(qtheta3, qbeta3);
    q04 = Quaternion_product(qtheta4, qbeta4);

    inv_q01 = Quaternion_inverse(q01);
    inv_q02 = Quaternion_inverse(q02);
    inv_q03 = Quaternion_inverse(q03);
    inv_q04 = Quaternion_inverse(q04);

    Op1 = Quaternion_product(Quaternion_product(q01, [0 L 0 0].'), inv_q01);
    Op1(1) = [];
    Op2 = Quaternion_product(Quaternion_product(q02, [0 L 0 0].'), inv_q02);
    Op2(1) = [];
    Op3 = Quaternion_product(Quaternion_product(q03, [0 L 0 0].'), inv_q03);
    Op3(1) = [];
    Op4 = Quaternion_product(Quaternion_product(q04, [0 L 0 0].'), inv_q04);
    Op4(1) = [];
    %% Find forces applied by all propellers thrusts on the body (in body frame)
    Tp1 = [0 0 0 kf*n1^2].';
    Tp2 = [0 0 0 kf*n2^2].';
    Tp3 = [0 0 0 kf*n3^2].';
    Tp4 = [0 0 0 kf*n4^2].';

    % P' = RPR'
    % P' = H(H(R, P), R')
    F1  = Quaternion_product(Quaternion_product(q1, Tp1), inv_q1);
    F2  = Quaternion_product(Quaternion_product(q2, Tp2), inv_q2);
    F3  = Quaternion_product(Quaternion_product(q3, Tp3), inv_q3);
    F4  = Quaternion_product(Quaternion_product(q4, Tp4), inv_q4);
    F1(1) = [];
    F2(1) = [];
    F3(1) = [];
    F4(1) = [];
    F = F1 + F2 + F3 + F4;

    %% Find torques applied by all propellers on the body (in body frame)
    Tauext1 = [0 0 0 km*n1^2].';
    Tauext2 = [0 0 0 -km*n2^2].';
    Tauext3 = [0 0 0 km*n3^2].';
    Tauext4 = [0 0 0 -km*n4^2].';
    % neglect rotor counter torque:
    % Tauext1 = [0 0 0 0].';
    % Tauext2 = [0 0 0 0].';
    % Tauext3 = [0 0 0 0].';
    % Tauext4 = [0 0 0 0].';
    
    M1  = Quaternion_product(Quaternion_product(q1, Tauext1), inv_q1);
    M2  = Quaternion_product(Quaternion_product(q2, Tauext2), inv_q2);
    M3  = Quaternion_product(Quaternion_product(q3, Tauext3), inv_q3);
    M4  = Quaternion_product(Quaternion_product(q4, Tauext4), inv_q4);
    M1(1) = [];
    M2(1) = [];
    M3(1) = [];
    M4(1) = [];
    M = M1 + M2 + M3 + M4 + cross(Op1, F1) + cross(Op2, F2) ...
        + cross(Op3, F3) + cross(Op4, F4);
    
    %% Find the force static matix 
    % matrix that link the force applied on the body to the vector:
    % kf*[cos(alpha1)n1^2; sin(alpha1)n1^2; cos(alpha2)n2^2; sin(alpha2)n2^2;
    %     cos(alpha3)n3^2; sin(alpha3)n3^2; cos(alpha4)n4^2; sin(alpha4)n4^2]
    [row columns] = size(n);
    A = [];
    Astatic = [];
    for ii = 1:3
        for jj = 1:columns
            Aij = F(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , ...
                {1,1,0}); Aij = subs(Aij,n , [0 0 0 0]);
            Ai2j = F(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)} , ...
                {1,1,pi/2}); Ai2j = subs(Ai2j,n , [0 0 0 0]);
            A = [A, Aij, Ai2j];
        end
         Astatic = [Astatic; A];
         A = [];
    end

    %% Find the torque static matix 
    % matrix that link the torque applied on the body to the vector:
    % kf*[cos(alpha1)n1^2; sin(alpha1)n1^2; cos(alpha2)n2^2; sin(alpha2)n2^2;
    %     cos(alpha3)n3^2; sin(alpha3)n3^2; cos(alpha4)n4^2; sin(alpha4)n4^2]
    A = [];
    A_M_stat = [];
    for ii = 1:3
        for jj = 1:columns
            Aij = M(ii); Aij = subs(Aij,{n(jj), kf, alpha(jj)} , {1,1,0});
            Aij = subs(Aij,n , [0 0 0 0]); Aij = subs(Aij,km , km/kf);
            Ai2j = M(ii); Ai2j = subs(Ai2j,{n(jj), kf, alpha(jj)} , {1,1,pi/2});
            Ai2j = subs(Ai2j,n , [0 0 0 0]); Ai2j = subs(Ai2j,km , km/kf);
            A = [A, Aij, Ai2j];
        end
         A_M_stat = [A_M_stat; A];
         A = [];
    end

    %% Test to verify solution
    
%     bRp1 = Quaternion_to_rotmx(q1);
%     bRp2 = Quaternion_to_rotmx(q2);
%     bRp3 = Quaternion_to_rotmx(q3);
%     bRp4 = Quaternion_to_rotmx(q4);
%     bRp1 = subs(bRp1, beta, beta0);
%     bRp1 = subs(bRp1, theta, theta0)
%     bRp2 = subs(bRp2, beta, beta0);
%     bRp2 = subs(bRp2, theta, theta0)
%     bRp3 = subs(bRp3, beta, beta0);
%     bRp3 = subs(bRp3, theta, theta0)
%     bRp4 = subs(bRp4, beta, beta0);
%     bRp4 = subs(bRp4, theta, theta0)
    
%     Op1 = subs(Op1, beta, beta0);
%     Op1 = subs(Op1, theta, theta0);
%     Op1 = subs(Op1, L, 1);
%     Op1 = round(Op1*10^5)/10^5
%     Op2 = subs(Op2, beta, beta0);
%     Op2 = subs(Op2, theta, tet
%     Op2 = round(Op2*10^5)/10^5
%     Op3 = subs(Op3, beta, beta0);
%     Op3 = subs(Op3, theta, theta0);
%     Op3 = subs(Op3, L, 1);
%     Op3 = round(Op3*10^5)/10^5
%     Op4 = subs(Op4, beta, beta0);
%     Op4 = subs(Op4, theta, theta0);
%     Op4 = subs(Op4, L, 1);
%     Op4 = round(Op4*10^5)/10^5
%     F = subs(F, beta, beta0);
%     F = subs(F, theta, theta0);
%     F = subs(F, alpha, alpha0);
%     Astatic1 = subs(Astatic, beta, beta0);
%     Astatic1 = subs(Astatic1, theta, theta0);
%     Astatic1 = double(round(Astatic1*10^5)/10^5)
%     
%     M = subs(M, beta, beta0);
%     M = subs(M, theta, theta0);
%     M = subs(M, alpha, alpha0);
%     A_M_stat = subs(A_M_stat, beta, beta0);
%     A_M_stat = subs(A_M_stat, theta, theta0)
end
% angle12 = rad2deg(double(round(atan2(norm(cross(Op1,Op2)), dot(Op1,Op2))*10^5)/10^5))
% angle23 = rad2deg(double(round(atan2(norm(cross(Op2,Op3)), dot(Op2,Op3))*10^5)/10^5))
% angle34 = rad2deg(double(round(atan2(norm(cross(Op3,Op4)), dot(Op3,Op4))*10^5)/10^5))
% angle41 = rad2deg(double(round(atan2(norm(cross(Op4,Op1)), dot(Op4,Op1))*10^5)/10^5))
% cross12 = cross(Op1,Op2)
% cross23 = cross(Op2,Op3)
% cross34 = cross(Op3,Op4)
% cross34 = cross(Op3,Op4)

