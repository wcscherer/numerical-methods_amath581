clear all; clc; close all;
% HW1 Initial Value Problems

% Consider the ODE dy/dt = -3ysin(t), y(0) = pi/sqrt(2)
% Exact solution is y(t) = (pi/sqrt(2)) * exp(3*cos(t)-3)

% 1a Solve the ODE numerically using the forward Euler method:

Ft=@(t) (pi/sqrt(2)).*exp(3.*cos(t)-3); %analytic solution
ft=@(y,t) -3.*y.*sin(t);                %dy/dt as f(t)

%Initialize the delt, t, y, E and meanE vectors
for i = 1:7
    delt(i) = 2^-(i+1);      % delta t step size for each simultaion
    t{i} = [0:2^-(i+1):5];   % vector of all t values for each simulation
    yE{i} = [pi/sqrt(2)];    % y vectors for each Euler simulation
    yH{i} = [pi/sqrt(2)];    % y vectors for each Heun simulation
    meanE(i) = 0;            % vector for mean error of each delt value
    meanH(i) = 0;
end

% forward euler and Heun method for all 7 simulations

for i = 1:7
    
    for j = 1:size(t{i},2)-1 % want to populate y up to the last time step
    
        yE{i}(j+1) = yE{i}(j) + delt(i).*ft(yE{i}(j),t{i}(j)); %forward Euler
        
        
        
        %constants for Heuns method for part 1b
        c1 = ft(yH{i}(end),t{i}(j));
        c2 = ft(yH{i}(end)+delt(i).*ft(yH{i}(end),t{i}(j)),t{i}(j)+delt(i));
        
        %implementing Heun's method using constants
        yH{i}(j+1) = yH{i}(j) + (delt(i)./2).*(c1+c2);
        
        
    end
    
    meanE(i) = mean(abs(Ft(t{i})-yE{i})); % finding mean error of all step errors for each dt
    
    meanH(i) = mean(abs(Ft(t{i})-yH{i})); % finding mean error of all step errors for each dt
end

%polyfit for mean error
p = polyfit(log10(delt), log10(meanE),1);

% %plot the error for each value of delta t
% figure(1)
% plot(log10(delt),log10(meanE),'x',log10(delt), p(1).*log10(delt)+p(2))
% title('Euler Method Error Analysis for dy/dt = -3ysin(t)')
% xlabel('log10 of Timestep')
% ylabel('log10 of Mean Error per Timestep')
% legend({'log10(Error)','Line Fit Slope = '+string(p(1))}, 'Location', 'northwest')

A1 = yE{7}'; A3 = p(1); % need to save values as variables so save works

save A1.dat  A1  -ascii % saving solution for for last euler scheme
save A2.dat  meanE -ascii % saving mean error for each euler scheme
save A3.dat  A3  -ascii % saving slope of poly fit of log error


% 1 b) Plot and save the solved ODE using Heun's method

%polyfit for mean error
u = polyfit(log(delt), log(meanH),1);

% %plot the error for each value of delta t
% figure(2)
% plot(log(delt),log(meanH),'x',log(delt), u(1).*log(delt)+u(2))
% title('Heun Method Error Analysis for dy/dt = -3ysin(t)')
% xlabel('log10 of Timestep')
% ylabel('log10 of Mean Error per Timestep')
% legend({'log10(Error)','Line Fit Slope = '+string(u(1))}, 'Location', 'northwest')

A4 = yH{7}.'; A6 = u(1); % need to save values as variables so save works

save A4.dat  A4  -ascii % saving solution for for last Heun scheme
save A5.dat  meanH -ascii % saving mean error for each Heun scheme
save A6.dat  A6  -ascii % saving slope of poly fit of log error

% problem 2

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
    
    t113mean(i) = mean(diff(soln113{i}.x));%average step size needed for each TOL of ode 113
    
end

slope_45 = polyfit(log10(t45mean),log10(TOL),1);
slope_23 = polyfit(log10(t23mean),log10(TOL),1);
slope_113 = polyfit(log10(t113mean),log10(TOL),1);

A8 = slope_45(1); A9 = slope_23(1); A10 = slope_113(1);

save A8.dat  A8  -ascii % saving ode45 truncation error
save A9.dat  A9  -ascii % saving ode23 truncation error
save A10.dat A10 -ascii % saving ode113 truncation error
