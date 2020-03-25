clear all; close all; clc;
%HW3 problem 1
xspan = [-1 1];
yspan = xspan;
npoints = 8;

%1a develop Dx and Dxx operator matrices with and without periodic boundary
%conditions A1 - A4

[A1 h] = Dx_1D(xspan,npoints,0); %Dx with no boundary coniditons 1D

A1 = full(A1);
save A1.dat  A1  -ascii

[A2 h] = Dx_1D(xspan,npoints,1); %Dx with periodic boundary conditions 1D

A2 = full(A2);
save A2.dat  A2  -ascii

[A3 h] = Dxx_1D(xspan,npoints,0); %Dxx with no boundary coniditons 1D

A3 = full(A3);
save A3.dat  A3  -ascii

[A4 h] = Dxx_1D(xspan,npoints,1); %Dxx with boundary coniditons 1D

A4 = full(A4);
save A4.dat  A4  -ascii

[A5, dx, dy] = Dx_2D(xspan, yspan, npoints, 0); %Dx without BCs 2D

A5 = full(A5);
save A5.dat  A5  -ascii

[A6, dx, dy] = Dx_2D(xspan, yspan, npoints, 1); %Dx with BCs 2D

A6 = full(A6);
save A6.dat  A6  -ascii

[A7, dx, dy] = Dy_2D(xspan, yspan, npoints, 0); %Dy without BCs 2D

A7 = full(A7);
save A7.dat  A7  -ascii

[A8, dx, dy] = Dy_2D(xspan, yspan, npoints, 1); %Dy with BCs 2D

A8 = full(A8);
save A8.dat  A8  -ascii

[A9, dx, dy] = Dxx_2D(xspan, yspan, npoints, 0); %Dxx without BCs 2D

A9 = full(A9);
save A9.dat  A9  -ascii

[A10, dx, dy] = Dxx_2D(xspan, yspan, npoints, 1); %Dxx with BCs 2D

A10 = full(A10);
save A10.dat  A10  -ascii

[A11, dx, dy] = Dyy_2D(xspan, yspan, npoints, 0); %Dyy without BCs 2D

A11 = full(A11);
save A11.dat  A11  -ascii

[A12, dx, dy] = Dyy_2D(xspan, yspan, npoints, 1); %Dyy with BCs 2D

A12 = full(A12);
save A12.dat  A12  -ascii

[A13, dx, dy] = Lxy_2D(xspan, yspan, npoints, 0); %Lxy without BCs 2D

A13 = full(A13);
save A13.dat  A13  -ascii

[A14, dx, dy] = Lxy_2D(xspan, yspan, npoints, 1); %Lxy with BCs 2D

A14 = full(A14);
save A14.dat  A14  -ascii

% Problem 2: Iterative solutions using Jacobi and Gauss-Seidel Methods

% initiate initial conditions
b10 = ones(10,1);

b30 = ones(30,1);

d10 = ones(10,1);

B10 = spdiags([-1.16*d10 d10 0.16*d10],-1:1,10,10);

d30 = ones(30,1);

B30 = spdiags([-1.16*d30 d30 0.16*d30],-1:1,30,30);

x10 = zeros(10,1);

x30 = zeros(30,1);

tol = 10^-8;
% Part a: Jacobi Method where M is the diagonal of B##

[xa10 A15 itra10] = Jacobi(x10,B10,b10,tol);

save A15.dat  A15  -ascii

[xa30 A16 itra30] = Jacobi(x30,B30,b30,tol);

save A16.dat  A16  -ascii

% Part b: Gauss-Seidel Method where M is the lower triangle of B##

[xb10 A17 itrb10] = Gauss_Seidel(x10,B10,b10,tol);

save A17.dat  A17  -ascii

[xb30 A18 itrb30] = Gauss_Seidel(x30,B30,b30,tol);

save A18.dat  A18  -ascii