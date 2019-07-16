function [Sign] = Sign(Arg)
%SIGN returns the sign and considers the 0 as positive
if Arg < 0
    Sign = -1;
else
    Sign =1;
end
end

