function [alphastar,wstar] = Mav_mirror_points(i, n,  D_unit2, d, alphastarHistoric, wstarHistoric, alphastar, wstar, opt_iterations)
%MAV_MIRROR_POINTS Summary of this function goes here
%   Detailed explanation goes here

%Look if the force, torque, hover space to see if it is symetric.
% if for better solutions in
%oppoit directions
dd = find(D_unit2(1,:) == - d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(2) = -Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(3) = -Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            alpha(4) = -Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(2) = Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(3) = -Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            alpha(4) = Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(2) = -Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(3) = Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            alpha(4) = -Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == -d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha = -alphastarHistoric(dd,:);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -alphastarHistoric(dd,1);
            alpha(2) = alphastarHistoric(dd,2);
            alpha(3) = -alphastarHistoric(dd,3);
            alpha(4) = alphastarHistoric(dd,4);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(2) = Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(3) = Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            alpha(4) = Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(1));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = alphastarHistoric(dd,1);
            alpha(2) = -alphastarHistoric(dd,2);
            alpha(3) = alphastarHistoric(dd,3);
            alpha(4) = -alphastarHistoric(dd,4);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, wstarHistoric(:,dd), L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(wstarHistoric(:,dd)*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = wstarHistoric(:,dd);
                    Hstar(i) = m*g/(kf*norm(wstarHistoric(:,dd))^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
%%%Opposit
dd = find(D_unit2(1,:) == - d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(2) = -Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(3) = -Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            alpha(4) = -Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(2) = Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(3) = -Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            alpha(4) = Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(2) = -Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(3) = Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            alpha(4) = -Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == -d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -alphastarHistoric(dd, 2);
            alpha(2) = -alphastarHistoric(dd, 1);
            alpha(3) = -alphastarHistoric(dd, 4);
            alpha(4) = -alphastarHistoric(dd, 3);
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -alphastarHistoric(dd, 2);
            alpha(2) = alphastarHistoric(dd, 1);
            alpha(3) = -alphastarHistoric(dd, 4);
            alpha(4) = alphastarHistoric(dd, 3);
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = Sign(alphastarHistoric(dd, 2))*(pi - abs(alphastarHistoric(dd, 2)));
            alpha(2) = Sign(alphastarHistoric(dd, 1))*(pi - abs(alphastarHistoric(dd, 1)));
            alpha(3) = Sign(alphastarHistoric(dd, 4))*(pi - abs(alphastarHistoric(dd, 4)));
            alpha(4) = Sign(alphastarHistoric(dd, 3))*(pi - abs(alphastarHistoric(dd, 3)));
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == -d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = alphastarHistoric(dd, 2);
            alpha(2) = -alphastarHistoric(dd, 1);
            alpha(3) = alphastarHistoric(dd, 4);
            alpha(4) = -alphastarHistoric(dd, 3);
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
dd = find(D_unit2(1,:) == d(2));
DD = D_unit2(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isequal(D_unit2(:,end), D_unit2(:,dd))
    if ~isempty(dd)
        if round(Hstar(i)*10^3)/10^3 < round(Heff(dd)*10^3)/10^3
            h1 = Hstar(:,i);
            alpha(1) = -alphastarHistoric(dd, 2);
            alpha(2) = -alphastarHistoric(dd, 1);
            alpha(3) = -alphastarHistoric(dd, 4);
            alpha(4) = -alphastarHistoric(dd, 3);
            n(1) = wstarHistoric(2,dd);
            n(2) = wstarHistoric(1,dd);
            n(3) = wstarHistoric(4,dd);
            n(4) = wstarHistoric(3,dd);
            % calculate angular and linear acceleration with this alphastar and wstar
            [m, ~, pdotdot, ~] = Quadcopter_tilted_arms_dynamic(kf, km, wRb, alpha, beta, theta, n, L, g, Mb, Mp, R, false);
            pdotdot = round(dec*pdotdot)/dec;
            Htest = m*pdotdot; % Force applied to the body with the propellers in this
            % if this solution does not break the constraint Mstar // d
            if isequal(round(Htest*10^2)/10^2,round(Fdes*10^2)/10^2) && i<opt_iterations
                test = [];
                [rows, ~] = size(alphastar);
                for j = 1:rows
                    if isequal(round(alpha*10^3)/10^3,round(alphastar(j,:)*10^3)/10^3) && isequal(round(n.'*10^3)/10^3, round(wstar(:, j)*10^3)/10^3)
                        test = [test true];
                    end
                end
                if isempty(test)
                    i = i+1;
                    alphastar(i,:) = alpha;
                    wstar(:,i) = n;
                    Hstar(i) = m*g/(kf*norm(n)^2);
                    i = i+1;
                    return;
                end
            end
        end
    end
end
end

