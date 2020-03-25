clear all; clc; close all;
Ft=@(t) (pi/sqrt(2)).*exp(3.*cos(t)-3); %analytic solution
ft=@(t,y) -3.*y.*sin(t);                %dy/dt as f(t,y)



for i = 1:7
    delt(i) = 2^-(i+1);      % delta t step size for each simultaion
    t{i} = [0];
    yH{i} = [pi/sqrt(2)];    % y vectors for each Heun simulation
    
    HE{i} = [0];             % vector for heun error
    meanH(i) = 0;
end


for i = 1
    
    while t{i} <= 5
        
        t{i}
        c1 = ft(t{i}, yH{i}(end));
        c2 = ft(t{i}+delt(i), yH{i}(end)+delt(i).*ft(t{i}, yH{i}(end)));
        
        yH{i}(end+1) = yH{i}(end) + (delt(i)./2).*(c1 + c2);
        
        t{i} = t{i} + delt(i);
        
    end
    
     %meanH(i) = mean(abs(Ft(t{i})-yH{i}));
    
    
end

%u = polyfit(log10(delt),log10(meanH),1);
