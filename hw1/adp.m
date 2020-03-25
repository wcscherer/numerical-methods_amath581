function dydt = adp(t,y,eps)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
dydt = [y(2); -eps*(y(1)^2-1)*y(2)-y(1)];
end

