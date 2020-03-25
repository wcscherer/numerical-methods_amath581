function [dxx_matrix, dx, dy] = Dxx_2D(xspan,yspan, np,val)
%Creates a 2D d^2/dx^2 operator matrix
%   This function creates a 2D d^2/dy^2 operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(1)-xspan(2))/np;
dy = abs(yspan(1)-yspan(2))/np;

if val < 1 % create an operator matrix without periodic boundary
    
    OD = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    D  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2
    T = (1/(dx*dx))*(OD+OD'+D); % generate sparse matrix
    
    I = eye(np);
    
    dxx_matrix = kron(I,T);
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    OD = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    D  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2
    
    
    B = zeros(np);
    
    B(np,1) = 1; B(1,np) = 1; % periodic boundary conditions
    
    T = (1/(dx*dx))*(OD+OD'+D+B); % generate sparse matrix
    
    I = eye(np);
    
    dxx_matrix = kron(I,T); % Dxx with periodic boundaries
    
    
end

end

