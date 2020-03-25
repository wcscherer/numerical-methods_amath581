function [dx_matrix,dx] = Dx_1D(xspan,np,val)
%Creates a 1D d/dx operator matrix
%   This function creates a 1D d/dx operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(2)-xspan(1))/np;


if val < 1 % create an operator matrix without periodic boundary
    
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    dx_matrix = (1/(2*dx))*(-D+D'); % generate sparse upper and lower matrix
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    D = sparse(2:np,1:np-1,ones(1,np-1),np,np); % offdiagonal below diag
    
    BCL = sparse(1,np,-1,np,np); % periodic BC at the left point
    
    BCR = sparse(np,1,1,np,np); % periodic BC at the right point
    
    dx_matrix = (1/(2*dx))*(-D+D'+ BCR + BCL);
    
    
end


end

