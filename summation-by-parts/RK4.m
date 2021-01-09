t0 = 0;
tN = 1;
sigma = 0.08;
N = 50;
[u t] = RungeK4 (t0, tN, sigma, N);
u(:)



function [u t] = RungeK4 (t0, tN, sigma, N)
    t = t0:1/N:tN;
    syms x;
    u = [exp(-0.5*((x-.5)/sigma)^2) NaN(1,length(t)-1)];
    h = (tN-t0)/N; % time step
        
    % initialize vars
    ua = 0;
    ub = 0;
    uc = 0;
    
    for i=1:length(t)-1
        ua = u(i) + 0.5*h*diff(-u(i));
        ub = u(i) + 0.5*h*diff(-ua);
        uc = u(i) + h*diff(-ub);
        u(i+1) = u(i) + (1/6)*h*( diff(-u(i) - 2*ua - 2*ub - uc) );
    end
end


    
    
    
    
    
    