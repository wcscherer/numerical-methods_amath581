clear all; clc; close all;
% HW1 Initial Value Problems

% Consider the ODE dy/dt = -3ysin(t), y(0) = pi/sqrt(2)
% Exact solution is y(t) = (pi/sqrt(2)) * exp(3*cos(t)-1)

% 1a Solve the ODE numerically using the forward Euler method:

Ft=@(t) (pi/sqrt(2)).*exp(3.*cos(t)-3); %analytic solution
ft=@(t,y) -3.*y.*sin(t);                %dy/dt as f(t)

for i = 1:7
    
    delt(i) = 2^-(1+i);
    t{i} = [0:delt(i):5];
    yE{i} = [pi/sqrt(2) zeros(1,size(t{i},2)-1)];
    EE{i} = [zeros(1,size(t{i},2))];
    meanE(i) = 0;
end


for i = 1
   
    for j = 1:size(t{i},2)
        j
        yE{i}(j) = yE{i}(j-1) + delt(i).*ft(t{i}(j-1),yE{i}(j-1));
        EE{i} = abs(Ft(t{i}(j))-yE{i}(j));
    end

    
    meanE(i) = mean(EE{i});
    
end

% %polyfit for mean error
% p = polyfit(log10(delt), log10(meanE),1);
% 
% %plot the error for each value of delta t
% figure(1)
% plot(log10(delt),log10(meanE),'x',log10(delt), p(1).*log10(delt)+p(2))
% title('Euler Method Error Analysis for dy/dt = -3ysin(t)')
% xlabel('log10 of Timestep')
% ylabel('log10 of Mean Error per Timestep')
% legend({'log10(Error)','Line Fit Slope = '+string(p(1))}, 'Location', 'northwest')
% 
% A1 = yE{7}; A3 = p(1); % need to save values as variables so save works
% 
% save A1.dat  A1  -ascii % saving solution for for last euler scheme
% save A2.dat  meanE -ascii % saving mean error for each euler scheme
% save A3.dat  A3  -ascii % saving slope of poly fit of log error