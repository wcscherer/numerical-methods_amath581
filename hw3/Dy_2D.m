function [dx_matrix, dx, dy] = Dy_2D(xspan,yspan, np,val)
%Creates a 2D d/dy operator matrix
%   This function creates a 2D d/dx operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(1)-xspan(2))/np;
dy = abs(yspan(1)-yspan(2))/np;

if val < 1 % create an operator matrix without periodic boundary
    
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    P = (1/(2*dy))*(-D+D'); % generate sparse upper and lower matrix
    
    I = eye(np);
    
    dx_matrix = kron(P,I);
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); % offdiagonal below diag
    
    B = zeros(np);
    
    B(np,1) = 1; B(1,np) = -1; % periodic boundary conditions
        
    P = (1/(2*dy))*(-D+D');
    
    I = eye(np);
    
    dx_matrix = kron(P,I)+ (1/(2*dy))*kron(B,I);
    
end
end

