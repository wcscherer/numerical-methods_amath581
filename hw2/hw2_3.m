close all; clear all; clc;

% problem 3 - optimization problem for quantum oscillator
% Pxx -[1*x^2 - e]P = 0 solve for 5 eigenvalue/functions

x0 = [-4:4]; % span for x from -4 to 4 
eig_start = [1, 2.532, 4.1, 6.1, 8]; % beginning eigen value guess
%y0(1) = [1, sqrt(16-eig_start(1))];
%Y0(2) = eig_start;
eig_row = zeros(5,1); % vector to store all 5 eigenvaues
Y0 = [eig_start(1),sqrt(16-eig_start(1))];
%fobj = @(Y) RHS_h23(Y); % locks in param values
% initializing values for phi solution of each mode
for i = 1:5
   
   phi{i} = [];
    
end
options = optimoptions('fminunc','Display','off','Diagnostics','off');

for n = 1:5
    
    y0 = [eig_start(n),sqrt(16-eig_start(n))];
    
    fobj = @(Y) RHS_h23(Y); % locks in param values
    [eig,error,a,b] = fminunc(fobj, y0,options); % run our optimization
    itr(n) = b.iterations;
    eig_row(n) = eig(1);
    y0 = [(sqrt(16-eig_start(n))),eig_start(n)];
    
end

A13 = abs(eig_row);
A14 = abs(itr);


