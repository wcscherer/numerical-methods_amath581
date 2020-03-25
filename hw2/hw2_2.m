clc; clear all; close all;

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

save A7.dat  A7  -ascii % saving solution for modes 1 - 5
save A8.dat  A8  -ascii 
save A9.dat  A9  -ascii 
save A10.dat  A10  -ascii
save A11.dat  A11  -ascii
save A12.dat  A12  -ascii % saving eigenvalue solutions


%plot(xspan,-phi{1},xspan,-phi{2},xspan,-phi{3},xspan,-phi{4},xspan,phi{5})
%legend('m1', 'm2', 'm3', 'm4', 'm5')
