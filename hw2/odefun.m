function dydt = odefun(x,y,dummy,eig)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (x^2 - eig)*y(1);


end