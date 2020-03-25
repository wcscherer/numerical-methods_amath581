function [error, t, y] = RHS_h23(Y0)
%RHS function for problem 2-3

xspan = [-4:0.1:4];
[t, y] = ode45(@(x,y) [y(2); (x^2 - Y0(1))*y(1)],xspan,[1;Y0(2)]);
 
 error = abs(y(end,2) + sqrt(16 - Y0(1))*y(end,1));
 
end

