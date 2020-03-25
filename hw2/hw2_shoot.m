function rhs = hw2_shoot(xspan,x,y,dummy,eig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rhs = [y(2),; (x.^2-eig).*y(1)]
end

