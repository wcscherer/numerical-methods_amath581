function rhs = phi2(y,eig)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
rhs = [y(2); (x.^2-eig).*y(1)];
end
