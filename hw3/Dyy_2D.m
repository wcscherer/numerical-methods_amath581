function [dyy_matrix, dx, dy] = Dyy_2D(xspan,yspan, np,val)
%Creates a 2D d^2/dy^2 operator matrix
%   This function creates a 2D d^2/dy^2 operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(1)-xspan(2))/np;
dy = abs(yspan(1)-yspan(2))/np;

if val < 1 % create an operator matrix without periodic boundary
    
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    C = sparse(1:np,1:np,-2*ones(1,np),np,np); % generate diagonal matrix
    
    P = (D+D'); % generate sparse upper and lower matrix
    
    I = eye(np);
    
    dyy_matrix = (1/dy^2)*(kron(C,I)+kron(P,I));
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    C = sparse(1:np,1:np,-2*ones(1,np),np,np); % generate diagonal matrix
    
    P = (D+D'); % generate sparse upper and lower matrix
    
    I = eye(np);
    
    B = zeros(np);
    
    B(np,1) = 1; B(1,np) = 1; % periodic boundary conditions
    
    dyy_matrix = (1/dy^2)*(kron(C,I)+kron(P,I)+kron(B,I));
            
    
end


end

