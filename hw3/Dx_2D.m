function [dx_matrix, dx, dy] = Dx_2D(xspan,yspan, np,val)
%Creates a 2D d/dx operator matrix
%   This function creates a 2D d/dx operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(1)-xspan(2))/np;
dy = abs(yspan(1)-yspan(2))/np;

if val < 1 % create an operator matrix without periodic boundary
    
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    P = (1/(2*dx))*(-D+D'); % generate sparse upper and lower matrix
    
    I = eye(np);
    
    dx_matrix = kron(I,P);
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); % offdiagonal below diag
    
    B = zeros(np);
    
    B(np,1) = 1; B(1,np) = -1; % periodic boundary conditions
        
    P = (1/(2*dx))*((-D+D')+B); % create ful matrix with BCs
    
    I = eye(np);
    
    dx_matrix = kron(I,P);
    
end
end

