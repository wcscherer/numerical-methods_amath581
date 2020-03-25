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

% plot(xspan,phi{1},col(1),xspan,phi{2},col(2),xspan,phi{3},col(3),xspan,phi{4},col(4),xspan,phi{5},col(5))
% legend('mode = '+string(1) + ' eig = '+string(eig_row(1)), 'mode = '+string(2) + ' eig = '+string(eig_row(2)),'mode = '+string(3) + ' eig = '+string(eig_row(3)), 'mode = '+string(4) + ' eig = '+string(eig_row(4)), 'mode = '+string(5) + ' eig = '+string(eig_row(5)))
% 
