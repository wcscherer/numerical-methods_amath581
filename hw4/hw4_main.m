% Homework 4 Due 11/28 at 1am
clc; clear all; close all;
% ut = k*uxx t = [0, 1], x = [0, 1], k = 0.1
% Use periodic boundary conditions u(t,0) = u(t,1) and intial heat
% distribution and u(0,x) = sin(2x)

% Problem 1 implement foreward Euler method Dt*ut = k*Dxx*u

k = 0.1; dt = 0.05; dx = 0.1;
tspan = [0:dt:1]; xspan = [0 1]; x0 = [0:dx:1];

u0=@(x) sin(2*pi*x); % 

[M1 range] = Dxx_1D(xspan,dx,0);
M1(1,10) = 1/dx^2; M1(11,2) = 1/dx^2;

U0 = u0(x0);
Ut = [U0'];

A1 = Fwd_Eul(M1,U0,tspan,k,dt);

A1 = A1(:,end); % save the last vector soution U(x,t=1)

save A1.dat A1 -ascii 

% Problem 2, reverse Euler method
%DtUn = kDxxUn+1

[M2 range] = Dxx_1D(xspan,dx,0);
M2(1,10) = 1/dx^2; M2(11,2) = 1/dx^2;

U0 = u0(x0);
Ut2 = [U0'];

A2 = Back_Eul(M2,U0,tspan,k,dt);

A2 = A2(:,end);
save A2.dat A2 -ascii 
% problem 2 part b dt = 0.01 and dx = 0.01
dt = 0.01; dx = 0.01;

tspan = [0:dt:1]; xspan = [0 1]; x0 = [0:dx:1];
N = size(x0,2);
[M3 range] = Dxx_1D(xspan,dx,0);
M3(1,N-1) = 1/dx^2; M3(N,2) = 1/dx^2;
U20 = u0(x0);

A3 = Back_Eul(M3,U20,tspan,k,dt);
A3 = A3(:,end);
save A3.dat A3 -ascii 

%Problem 3: Crank Nicholson Method DtUn = k/2(DxxUn + DxxUn+1)

% part a, dx = 0.1 and dt = 0.05; A4 should be 11x1
dt = 0.05; dx = 0.1;
tspan = [0:dt:1]; xspan = [0 1]; x0 = [0:dx:1];
N = size(x0,2);
[M4 range] = Dxx_1D(xspan,dx,0);
M4(1,N-1) = 1/dx^2; M4(N,2) = 1/dx^2;

U0 = u0(x0);
Ut = [U0'];

A4 = Crk_Nln(M4,U0,tspan,k,dt);
A4 = A4(:,end);
save A4.dat A4 -ascii 

% part b, dx = 0.01 and dt = 0.01; A4 should be 11x1
dt = 0.01; dx = 0.01;
tspan = [0:dt:1]; xspan = [0 1]; x0 = [0:dx:1];
N = size(x0,2);
[M5 range] = Dxx_1D(xspan,dx,0);
M5(1,N-1) = 1/dx^2; M5(N,2) = 1/dx^2;

U0 = u0(x0);
Ut = [U0'];

A5 = Crk_Nln(M5,U0,tspan,k,dt);
A5 = A5(:,end);
save A5.dat A5 -ascii 