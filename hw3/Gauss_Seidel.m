function [xout,error,itr] = Gauss_Seidel(xinp,A, b, tol)
%Function implements the Gauss-Seidel method on matrix A
%   Implementing the Gauss-Seidel method on A and returns the vector xout
%   that solves the statement Axk = b, the error of the solution vector
%   based on the tolerance tol

itr = 0;
xk = xinp;
er = norm(b-A*xinp);
error = [er];
while er > tol & itr < 999
   
    % xk+1 = xk + M\r
    
   M = tril(full(A));% GS uses the lower triangle of A as the preconditioner
   
   r = (b - A*xk);
    
   xk =  xk + M\r;
    
   er = norm(b - A*xk);                                           
   
   itr = itr + 1;
   
   error = [error, er];
end

xout = xk;

itr = size(error,2);

itr;
end

