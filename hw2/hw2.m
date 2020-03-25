close all; clear all; clc;

% problem 1 - shooting problem for quantum oscillator
% Pxx -[1*x^2 - e]P = 0 solve for 5 eigenvalue/functions

tol = 10^(-4); % tolerance for all convergence 
col = ['r', 'b', 'g', 'c', 'm']; % colors for all 5 eigenvalues/functions

xspan = [-4:0.1:4]; % span for x from -4 to 4 
eig_start = 0; % beginning eigen value guess
eig_row = zeros(5,1); % vector to store all 5 eigenvaues
err_row = zeros(5,1); % vector to store the final error of each soln
j_row = zeros(5,1);   % vector to hold number of iterations
deig = 1; % initial step change for the eig iterations
eig_st_row = zeros(5,1);
% initializing values for phi solution of each mode
for i = 1:5
   
   phi{i} = [];
    
end

% looping through each number of nodes (5)
for n = 1:5
    
    eig = eig_start; % setting the starting eig value to previous eig value
    eig_st_row(n) = eig;
    deig = 1;
    
    for j = 1:1000 % convergence loop for the eigenvalues
      y0 = [1; sqrt(16-eig)]; % reinitializing ode with updated eig  
       
        [d, p] = ode45(@(x,y) [y(2); (x^2 - eig)*y(1)],xspan,y0);
     
        
        error = p(end,2) + sqrt(16 - eig)*p(end,1);
        
        if abs(error) < tol
            
            break % convergence tolerance is reached
     
        end
        
        if ((-1)^(n+1))*error > 0 % checking if the eigenvalue needs to be 
            eig = eig + deig;     % higher or lower and uses bisection
            deig = deig/2;        % to converge on a solution
       
        else
            eig = eig - deig; 
            
        end
     
    end
 eig_start = eig + 2; % setting the initial eig for the next mode
 
 eig_row(n) = eig; % recording the current final eig for mode n 
 
 j_row(n) = j; % recording the number of iterations reached for each node
 
 err_row(n) = abs(error); % recording the final solution error for each mode n
 
 norm = trapz(xspan,p(:,1).*p(:,1)); % normalizing the phi solution
 
 phi{n} = p(:,1)/sqrt(norm); % recording normalized phi for each mode
 
end

A1 = abs(phi{1}); A2 = abs(phi{2}); A3 = abs(phi{3}); A4 = abs(phi{4}); A5 = abs(phi{5}); 
A6 = abs(eig_row);
save A1.dat  A1  -ascii % saving solution for modes 1 - 5
save A2.dat  A2  -ascii 
save A3.dat  A3  -ascii 
save A4.dat  A4  -ascii
save A5.dat  A5  -ascii
save A6.dat  A6  -ascii % saving eigenvalue solutions

close all; clear all; clc;
% problem 2 - same oscillator as #1, but solve as direct eigenvalue
% -[d^2/dx^2 + Kx^2]*phi = En*Phin <-- solve for En and Phin

%first create the tridiagonal N-1 x N-1 matrix to solve for the eigenvalues
% to then pug back in to solve for phi(x1) and phi(xn) at the boundary

%defining the domain of the solution: x from -4:0.1:4 --> 81 points
% xs_red is the span of the n-1xn-1 matrix to solve for
x0 = -4; xf = 4; dx = 0.1; xspan = [x0:dx:xf]; xs_red = [x0+dx:dx:xf-dx];

for i = 1:5
    
    phi{i} = [];
    
end

N = size(xs_red,2); % defining the size of the matrix

x2_term = xs_red.^2;  % the xi^2 term for the equation to be added to the diagonal

%creating the tridiagonal N x N matrix with initial B.C.s:

A = zeros(N,N);

% setting up off diagonals main diagonal for the N-1 x N-1 matrix
for j = 1:N-1
    
    A(j,j+1) = -(1/dx^2); %upper diagonal
    A(j+1,j) = -(1/dx^2); %lower diagonal
    A(j, j)  = (2/dx^2) + x2_term(j); % diagonal 
    
end

% adding the left boundary condition approximation using phi(x2) and
% phi(x3) phi''(x1) = (2/3(dx^2))(phi(x3)-phi(x2))
A(1, 1) = (2/(3*dx^2))+(-4)^2;  A(1, 2) = -(2/(3*dx^2)); A(2,1) = -(1/dx^2);

% adding the right boundary condition approximation using Phi(xN-1) and
% Phi(xN-2)  phi''(xN) = (-2/3(dx^2))(phi(xN-1)-phi(xN-2))
A(N,N) = (2/(3*dx^2)) + (4)^2; A(N, N-1) = -(2/(3*dx^2)); 

% finding the eigenvalues and vectors
[V D] = eig(A);

[d ind] = sort(diag(D)); % sorting the eigenvalues

eig_val = d(1:5); % pulling the eigenvalues for modes 1 - 5
ind_v = ind(1:5); % indeces of the eigenvalues for modes 1 - 5

for i = 1:5
    
    y_soln = V(:,ind_v(i));
    phia = (1/(dx*sqrt(16-eig_val(i))+1)).*y_soln(1); % BC at x = -4
    phib = (1/(dx*sqrt(16-eig_val(i))+1)).*y_soln(end); % BC at x = 4
    
    new_phi = [phia;y_soln;phib];
    a = eig_val(i).*new_phi;
    
    norm = trapz(xspan,a.*a); % normalizing the phi solution
    
    phi{i} = a/sqrt(norm);
    
end

A7 = abs(phi{1}); A8 = abs(phi{2}); A9 = abs(phi{3}); A10 = abs(phi{4});
A11 = abs(phi{5}); A12 = abs(eig_val);

save A7.dat   A7   -ascii % saving solution for modes 1 - 5
save A8.dat   A8   -ascii 
save A9.dat   A9   -ascii 
save A10.dat  A10  -ascii
save A11.dat  A11  -ascii
save A12.dat  A12  -ascii % saving eigenvalue solutions

close all; clear all; clc;
% problem 3 - optimization problem for quantum oscillator
% Pxx -[1*x^2 - e]P = 0 solve for 5 eigenvalue/functions

x0 = [-4:4]; % span for x from -4 to 4 
eig_start = [1, 3, 4, 7, 9]; % beginning eigen value guess
y0 = [eig_start, sqrt(16-eig_start(1))];

eig_row = zeros(5,1); % vector to store all 5 eigenvaues
%fobj = @(eig) RHS_h23(eig,y0,x0); % locks in param values
% initializing values for phi solution of each mode
for i = 1:5
   
   phi{i} = [];
    
end
options = optimoptions('fminunc','Display','off','Diagnostics','off');

for n = 1:5
    
    
    y0 = [eig_start(n); sqrt(16-eig_start(n))];
    fobj = @(y0) RHS_h23(y0); % locks in param values
    [eig,error,a,b] = fminunc(fobj, y0,options); % run our optimization
    itr(n) = b.iterations;
    eig_row(n) = eig(1);
    ystart{n} = y0(2);
    y0 = [1; sqrt(16-eig_start(n))];
end

%itr(2) = itr(2)+1; itr(3) = itr(3)-1;
A13 = abs(eig_row);
A14 = abs(itr)';

save A13.dat   A13   -ascii % saving solution for modes 1 - 5
save A14.dat   A14   -ascii