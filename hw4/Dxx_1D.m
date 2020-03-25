function [dxx_matrix, span] = Dxx_1D(xspan,dx,val)
%Creates a 1D d^2/dx^2 operator matrix
%   This function creates a 1D d^2/dx^2 operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions

span = [xspan(1):dx:xspan(2)];
np = size(span,2);


if val < 1 % create an operator matrix without periodic boundary
    
    OD = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    D  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2
    dxx_matrix = (1/(dx*dx))*(OD+OD'+D); % generate sparse matrix
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    OD = sparse(2:np,1:np-1,ones(1,np-1),np,np); % offdiagonal below diag
    D  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2
    BCL = sparse(1,np,1,np,np); % periodic BC at the left point
    
    BCR = sparse(np,1,1,np,np); % periodic BC at the right point
    
    dxx_matrix = (1/(dx*dx))*(OD+OD'+ D + BCR + BCL);
    
    
end
end

% periodic BC at the left point