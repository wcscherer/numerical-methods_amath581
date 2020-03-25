function [Lxy_matrix, dx, dy] = Lxy_2D(xspan,yspan, np,val)
%Creates a 2D d^2/dx^2 + d^2/dy^2 laplacian operator matrix
%   This function creates a 2D d^2/dx^2 + d^2/dy^2 laplacian operator with or without
%   periodic boundary conditions over a number of grid points np
%   of a span xspan and yspan

%If val = 0: non-periodic boundary conditions
%If val = 1:  periodic boundary conditions


dx = abs(xspan(1)-xspan(2))/np;
dy = abs(yspan(1)-yspan(2))/np;

if val < 1 % create an operator matrix without periodic boundary
    I = eye(np);
    % calculating Dxx matrix
    ODx = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    Dx  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2 Dxx
    T = (1/(dx*dx))*(ODx+ODx'+Dx); % generate sparse matrix
    
    %calculating Dyy matrix
    ODy = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    Dy = (1/dy^2)*sparse(1:np,1:np,-2*ones(1,np),np,np); % generate diagonal Dyy matrix
    
    S = (1/dy^2)*(ODy+ODy'); % generate sparse upper and lower matrix
    
    
    Lxy_matrix = kron(I,T) + kron(S,I) + kron(I,Dy); % Lxy matrix
    
end



if val >= 1 % create an operator matrix with periodic boundary
     
    I = eye(np);
    % calculating Dxx matrix
    ODx = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    Dx  = sparse(1:np,1:np,-2*(ones(1,np)),np,np); % diagonal of -2 Dxx
    
    Bx = zeros(np);
    
    Bx(np,1) = 1; Bx(1,np) = 1; % periodic boundary conditions
    
    T = (1/(dx*dx))*(ODx+ODx'+Dx+Bx); % generate sparse matrix
    
    %calculating Dyy matrix
    ODy = sparse(2:np,1:np-1,ones(1,np-1),np,np); %generate sparse lower diag
    
    Dy = (1/dy^2)*sparse(1:np,1:np,-2*ones(1,np),np,np); % generate diagonal Dyy matrix
    
    S = (1/dy^2)*(ODy+ODy'); % generate sparse upper and lower matrix
    
    By = zeros(np);
    
    By(np,1) = 1; By(1,np) = 1; % periodic boundary conditions for y
    
    By = (1/dy^2)*By;
    
    Lxy_matrix = kron(I,T) + kron(S,I) + kron(I,Dy)+kron(By,I); % Lxy matrix
    
    
end
end

