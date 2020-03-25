clear all; clc; close all;
%HW 5 Spectral Methods

%Problem 1: solve the BVP Uxx + 4Ux + exp(x)U = sin(8x)
%Over X = [-1 1] and U(+-1) = 0; save U(0) as A1.dat

%Set up the Chedbychev Differential Operator Dx using the cheb function 
%from trefethen - Solve Lu = f

% L = Dxx + 4*Dx + exp(x)*I;  f = sin(8x)

N = 40;                    % Number of points modeled
[Dx,xspan] = cheb(N);      % Create the chebychev derivative matrix and xspan
Dxx = Dx^2;
Dx  = Dx(2:N,2:N);
Dxx = Dxx(2:N,2:N);        % BCs are U(+-1) = 0
  
% Solve system of equations
L = (Dxx + 4*Dx + diag(exp(xspan(2:N))));  % Matrix L operating on U
f = sin(8*xspan(2:N));     % right hand side        
U = L\f;                   % Solve for U over i = 2 to i = N
U = [0;U;0];               % include boundaries (0)

A1 = U(21);   % save solution at x = 0 corresponding to location 21
save A1.dat A1 -ascii


% Problem 2: modify program 14 to use Newtons Method instead using
% chebyshev matrices to solve the non linear ODE: Uxx = exp(u)

% When using Newton's method, solve Dxx*U - exp(U) = 0
% Un+1 = Un - J\F(U) where J is the Jacobian of F

N = 16;
[Dx,x] = cheb(N); Dxx = Dx^2; Dxx = Dxx(2:N,2:N);
U = zeros(N-1,1);
itr = 0;
err = [1];
  while err > 1e-15
    F = Dxx*U - exp(U);          % F vector
  	J = Dxx - diag(exp(U));      % Jacobian matrix: d/du(F)
    
    Unew = U - J\F;              % Newton's method update of U
    err(end+1) = norm(Unew - U,inf);    % Calculate the euclidian norm 
    itr = itr + 1;
    U = Unew;
  end
  
U = [0;U;0];                    % Save final vector solution of U with boundaries

A2 = err(2); A3 = err(3);       % Save the first and second iteration step errors

save A2.dat A2 -ascii
save A3.dat A3 -ascii

% Problem 3: solve the non-linear BVP Ut = Dxx*U + exp(U) over x = [-1 1]
% with boundaries u(+-1,t) = 0 and U(x,0) = 0

% determining 
N = 40;                    % Number of points modeled
[Dx,xspan] = cheb(N);      % Create the chebychev derivative matrix and xspan
Dxx = Dx^2;
Dx  = Dx(2:N,2:N);
Dxx = Dxx(2:N,2:N);        % BCs are U(+-1) = 0
tspan = [0:10^-3:3.55];    % Need fine tspan 
U0 = zeros(N-1,1);

[t U] = ode23s(@(t,u,x) Dxx*u + exp(u),tspan,U0); % solve the PDE
A4 = U(3501,20); % Should be U(0,3.5)

U0T = U(:,20); %U(0,t)
A5 = interp1(U0T,tspan,5); % interpolate to find where t5 where U(0,t5) = 5;


save A4.dat A4 -ascii
save A5.dat A5 -ascii
