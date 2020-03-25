clc; clear all; close all;
% question 2 of hw1
% solve y'' + eps(y^2 - 1)y' + y = 0

%re-write to first order y'1 = y2, y'2 = -eps(y1^2-1)y2-y1

% problem 2a
t2a = [0:0.5:32];
y0a = [sqrt(3);1];

% eps = 0.1
[t45_01,y45_01] = ode45(@(t,y) [y(2); -0.1*(y(1)^2-1)*y(2)-y(1)],t2a,y0a);


% eps = 1
[t45_1,y45_1] = ode45(@(t,y) [y(2); -1*(y(1)^2-1)*y(2)-y(1)],t2a,y0a);


% eps = 20
[t45_20,y45_20] = ode45(@(t,y) [y(2); -20*(y(1)^2-1)*y(2)-y(1)],t2a,y0a);

A7 = [y45_01(:,1), y45_1(:,1), y45_20(:, 1)];
save A7.dat  A7  -ascii % saving all three y1 solutions for all epsilons


% problem 2b

% initial conditions
y0b = [2,pi^2];
t2span = [0:32];

% use ode45 and eps = 1 with a tolerance to solve van der Pol eqn for
% various tolerances from 10^-4 to 10^-10

for i = 4:1:10
    
    TOL(i-3) = 10^-i;
    options(i-3) = odeset('AbsTol',TOL(i-3),'RelTol',TOL(i-3));
    soln45{i-3} = [];
    soln23{i-3} = [];
    soln113{i-3} = [];
    t45mean(i-3) = 0;
    t23mean(i-3) = 0;
    t113mean(i-3) = 0;
end

% solve the oscillayor using ode45, ode23, and ode113
for i = 1:size(TOL,2)
    
    % eps = 1
    soln45{i} = ode45(@(t,y) [y(2); -1*(y(1)^2-1)*y(2)-y(1)],t2span,y0b,options(i));

    soln23{i} = ode23(@(t,y) [y(2); -1*(y(1)^2-1)*y(2)-y(1)],t2span,y0b,options(i));
    
    soln113{i} = ode113(@(t,y) [y(2); -1*(y(1)^2-1)*y(2)-y(1)],t2span,y0b,options(i));

end

% calculate the average step size needed for each ode solver for a given
% tolerance tmean = mean(diff(T))

for i = 1:size(TOL,2)
    
    t45mean(i) = mean(diff(soln45{i}.x)); %average step size needed for each TOL for ode45
    
    t23mean(i) = mean(diff(soln23{i}.x)); %average steo size needed for each TOL for ode23
    
    t113mean(i) = mean(diff(soln23{i}.x));%average step size needed for each TOL of ode 113
    
end

slope_45 = polyfit(log10(t45mean),log10(TOL),1);
slope_23 = polyfit(log10(t23mean),log10(TOL),1);
slope_113 = polyfit(log10(t113mean),log10(TOL),1);

A8 = slope_45(1); A9 = slope_23(1); A10 = slope_113(1);

save A8.dat  A8  -ascii % saving ode45 truncation error
save A9.dat  A9  -ascii % saving ode23 truncation error
save A10.dat A10 -ascii % saving ode113 truncation error