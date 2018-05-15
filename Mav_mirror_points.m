function [alphai,wi] = Mav_mirror_points(n, number_of_directions, D_unit, d, alphastarHistoric, wstarHistoric, alphastar, wstar, to_check, checker, dec)
%MAV_MIRROR_POINTS Summary of this function goes here
%   Detailed explanation goes here
w = zeros(n,1);
alpha = zeros(n,1);
wi = [];
alphai = [];

%Look if the force, torque, hover space to see if it is symetric.
% if for better solutions in
%oppoit directions
dd = find(D_unit(1,:) == - d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                alpha(ll) = -Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
                else
                    alpha(ll) = -Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
                end
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
                else
                    alpha(ll) = Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
                end
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == -d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            alpha = -alphastarHistoric(:,dd);
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = alphastarHistoric(ll, dd);
                else
                    alpha(ll) = -alphastarHistoric(ll, dd);
                end
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                alpha(ll) = Sign(alphastarHistoric(ll, dd))*(pi - abs(alphastarHistoric(ll, dd)));
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -alphastarHistoric(ll, dd);
                else
                    alpha(ll) = alphastarHistoric(ll, dd);
                end
            end
            w = wstarHistoric(:,dd);
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
%%%Opposit
dd = find(D_unit(1,:) == - d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -Sign(alphastarHistoric(ll-1, dd))*(pi- abs(alphastarHistoric(ll-1, dd)));
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = -Sign(alphastarHistoric(ll+1, dd))*(pi- abs(alphastarHistoric(ll+1, dd)));
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = Sign(alphastarHistoric(ll-1, dd))*(pi- abs(alphastarHistoric(ll-1, dd)));
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = -Sign(alphastarHistoric(ll+1, dd))*(pi- abs(alphastarHistoric(ll+1, dd)));
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -Sign(alphastarHistoric(ll-1, dd))*(pi- abs(alphastarHistoric(ll-1, dd)));
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = Sign(alphastarHistoric(ll+1, dd))*(pi- abs(alphastarHistoric(ll+1, dd)));
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == -d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = -alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = -alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = Sign(alphastarHistoric(ll-1, dd))*(pi -abs(alphastarHistoric(ll-1, dd)));
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = Sign(alphastarHistoric(ll+1, dd))*(pi -abs(alphastarHistoric(ll+1, dd)));
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*dec/100)/dec/100 < round(checker(dd)*dec/100)/dec/100
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll < n
                        alpha(ll) = alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*dec/100)/dec/100,round(alphastar(:,j)*dec/100)/dec/100) && isequal(round(w*dec/100)/dec/100, round(wstar(:, j)*dec/100)/dec/100)
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end

%% Methode 2
dd = find(D_unit(1,:) == - d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                alpha(ll) = -Sign(alphastarHistoric(ll, dd))  +alphastarHistoric(ll, dd);
                w(ll) = wstarHistoric(ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = Sign(alphastarHistoric(ll-n, dd))  -alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = Sign(alphastarHistoric(ll, dd))  -alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            for ll = 2:2:n
                alpha(ll) = Sign(alphastarHistoric(ll, dd)) -alphastarHistoric(ll, dd);
                w(ll) = wstarHistoric(ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = Sign(alphastarHistoric(ll-n, dd))  -alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = Sign(alphastarHistoric(ll, dd))  -alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            for ll = 1:2:n
                alpha(ll) = Sign(alphastarHistoric(ll, dd)) -alphastarHistoric(ll, dd);
                w(ll) = wstarHistoric(ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == -d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = -alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = -alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = -alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = -alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            for ll = 1:2:n
                alpha(ll) = -alphastarHistoric(ll, dd);
                w(ll) = wstarHistoric(ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = -Sign(alphastarHistoric(ll-n, dd))*pi + alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = -Sign(alphastarHistoric(ll, dd))*pi + alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(1));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(2));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = -alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = -alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            for ll = 2:2:n
                alpha(ll) = -alphastarHistoric(ll, dd);
                w(ll) = wstarHistoric(ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
%%% Sides symetries
dd = find(D_unit(1,:) == - d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = Sign(alphastarHistoric(ll-1, dd))*pi - alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = Sign(alphastarHistoric(ll+1, dd))*pi  - alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -Sign(alphastarHistoric(ll-1, dd))*pi  + alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = Sign(alphastarHistoric(ll+1, dd))*pi  - alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = Sign(alphastarHistoric(ll-1, dd))*pi  - alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = -Sign(alphastarHistoric(ll+1, dd))*pi  + alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == -d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 3:n+2
                if ll > n
                    alpha(ll-2) = alphastarHistoric(ll-n, dd);
                    w(ll-2) = wstarHistoric(ll-n,dd);
                else
                    alpha(ll-2) = alphastarHistoric(ll, dd);
                    w(ll-2) = wstarHistoric(ll,dd);
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == - d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == - d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = -Sign(alphastarHistoric(ll-1, dd))*pi  + alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = -Sign(alphastarHistoric(ll+1, dd))*pi  + alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == -d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                if mod(ll,2) == 0
                    alpha(ll) = alphastarHistoric(ll-1, dd);
                    w(ll) = wstarHistoric(ll-1,dd);
                else
                    if ll ~=n
                        alpha(ll) = -alphastarHistoric(ll+1, dd);
                        w(ll) = wstarHistoric(ll+1,dd);
                    end
                end
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
dd = find(D_unit(1,:) == d(2));
DD = D_unit(:,dd);
dd2 = find(DD(2,:) == d(1));
DD = DD(:,dd2);
dd3  = DD(3,:) == d(3);
dd = dd(dd2(dd3));
if ~isempty(dd)
    if dd < number_of_directions
        if round(to_check*(dec/100))/(dec/100) < round(checker(dd)*(dec/100))/(dec/100)
            alpha = zeros(n,1);
            w = zeros(n,1);
            for ll = 1:n
                alpha(ll) = alphastarHistoric((n+1)-ll, dd);
                w(ll) = wstarHistoric((n+1)-ll,dd);
            end
            test = [];
            [~, columns] = size(alphastar);
            for j = 1:columns
                if isequal(round(alpha*(dec/100))/(dec/100),round(alphastar(:,j)*(dec/100))/(dec/100)) && isequal(round(w*(dec/100))/(dec/100), round(wstar(:, j)*(dec/100))/(dec/100))
                    test = [test true];
                end
            end
            if isempty(test)
                alphai = alpha;
                wi = w;
                return;
            end
        end
    end
end
wi = [];
alphai = [];
end


