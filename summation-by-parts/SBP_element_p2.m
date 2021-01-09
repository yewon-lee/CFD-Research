%% SBP using element-type operators, p=2 (Q1.6) ops 1.88-1.93

% Errorless to my knolwedge.
% Note: for element type ops, p = (# of nodes) - 1

% close all; clear; clc
% 
% a = 1;
% tN = 5;
% t0 = 0;
% x0 = 0;
% xN = 1;
% sigma0 = 0;
% sigma1 = 1;
% elements1 = 5;
% elements2 = 5;
% elements3 = 60;
% elements4 = 50;
% [u1 x1 err1 uHu1 oneHu1 t1 dx1] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements1);
% [u2 x2 err2 uHu2 oneHu2 t2 dx2] = SBP_solver (a, tN, t0, x0, xN, sigma1, elements2);
% % [u3 x3 err3 uHu3 oneHu3 t3 dx3] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements3);
% % [u4 x4 err4 uHu4 oneHu4 t4 dx4] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements4);
% 
% % %plot
% % plot(x1,u1);
% % xlabel('x');
% % ylabel('u');
% % hold on
% % xzero = 0:1/50:1;
% % xzero = xzero';
% % uzero = NaN(length(xzero),1);
% % for i=1:length(uzero)
% %     uzero(i) = exp(-0.5*((xzero(i)-0.5)/0.08)^2);
% % end
% % % hold on
% % % plot(x2,u2)
% % %plot(x2,u2)
% % %plot(x3,u3)
% % %plot(x4,u4)
% % plot(xzero, uzero);
% % legend ('t=1','t=5', 'exact soln');
% % title('u(x,5)');
% 
% % % error - loglog plot
% % deltax = [dx1 dx2 dx3 dx4]';
% % errvals = [err1 err2 err3 err4]';
% % loglog(deltax, errvals, 'x')
% % xlabel('delta x');
% % ylabel('error');
% % hold on
% % 
% % p = polyfit(log(deltax),log(errvals),1)
% % err_ = polyval(p,log(deltax));
% % loglog(deltax,exp(err_));
% % [m,b] = polyfit(log(deltax),log(errvals),1);
% 
% 
% %uTranspose * H * u
% plot(t1,uHu1);
% hold on
% plot(t2, uHu2)
% % title('symmetric energy uTranspose*H*u vs. t');
% legend('upwind','symm');
% 
% % % onesTranspose * H * u
% % %close all; clear; clc
% % plot(t1,oneHu1);
% % hold on
% % plot(t2,oneHu2);
% % % plot(t3,oneHu3);
% % % plot(t4,oneHu4);
% % title('onesTranspose*H*u vs. t');
% % legend('upwind','symmetric');

% ----------------------function---------------------

% Returns u(x) at t=tN
%
% Returns rate of convergence with grid refinement of the H-norm of the
% solution error
%
% Returns uTranspose * H * u as a function of time (= energy)
%
% Returns 1Transpose * H * u as a function of time (= conservation for e.g. of mass)
%
% Notes: sigma = 1 -> symmetric
%        sigma = 0 -> upwind
%      
%        elements is the number of elements
%
%

function[u x err uHu oneHu t dx] = SBP_element_p2 (a, tN, t0, x0, xN, sigma, elements)

    % x
    xi = [-sqrt(15)/5   0   sqrt(15)/5]; 
    x = NaN(3,elements);
    xmax = NaN(1,elements);
    xmin = NaN(1,elements);
    dElems = (xN - x0)/elements;
    for i = 1:elements
        xmin(i) = x0 + dElems*(i-1);
    end
    for i = 1:elements
        if i == elements
            xmax(i) = xN;
        else
            xmax(i) = xmin(i+1);
        end
    end
    for i = 1:elements
        x(:,i) = (xmax(i)-xmin(i)) / (1-(-1)) * (xi+1) +  xmin(i);
    end
    xElems = x;
    x = reshape(xElems,[3*elements,1]);
      
    % u init
    uElems = NaN(3, elements);
    for i = 1:3
        for j = 1:elements
            uElems(i,j) = exp(-0.5*( (xElems(i,j) - 0.5) / 0.08 )^2);
        end
    end
    
    % Dz matrix
    Dz = [-sqrt(15)/2   2*sqrt(15)/3   -sqrt(15)/6;
          -sqrt(15)/6   0              sqrt(15)/6;
          sqrt(15)/6   -2*sqrt(15)/3   sqrt(15)/2];
    J = ( (xN-x0) / elements ) / 2;
    D = Dz/J;
        
    % H matrix
    v = [5/9   8/9   5/9];
    Hz = zeros(3) + diag(v,0);
    H = J*Hz;
    
    % Q matrix
    Q = H*D;
           
    % elementary vectors 
    NL = 3;
    tL = [(sqrt(15)+5)/6   -2/3   (-sqrt(15)+5)/6]';
    tR = [(-sqrt(15)+5)/6   -2/3   (sqrt(15)+5)/6]';
    
    % case upwind    
    if sigma == 0
        
        % A matrix
        A = -a*(D + inv(H)*tL*tL');
            
        % B matrix
        B = inv(H)*a*tL*tR';     
        
        % C matrix
        C = zeros(NL);
        
        
    % case symmetric  
    elseif sigma == 1
        
        % A matrix
        A = (-D*a - 0.5*inv(H)*a*tL*tL' + 0.5*inv(H)*a*tR*tR');
                    
        % B matrix
        B = 0.5*inv(H)*a*tL*tR';
        
        % C matrix
        C = -0.5*inv(H)*a*tR*tL';

    end
        
    % time marching
    dx = dElems/4;
    dt = 0.1*dx;
    N = floor((tN-t0)/dt);
    dt = (tN-t0)/N;
    uhist = NaN(NL*elements, N);
    ua = NaN(NL,elements); 
    ub = NaN(NL,elements); 
    uc = NaN(NL,elements); 
    ufin = NaN(NL,elements); 
    for j = 1:N
        for i = 1:elements 
            if i == 1
                R = A*uElems(:,i) + B*uElems(:,elements)+ C*uElems(:,i+1);
                ua(:,i) = uElems(:,i) + 0.5 * dt * (A*uElems(:,i) + B*uElems(:,elements) + C*uElems(:,i+1));
            elseif i == elements
                R = A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,1);
                ua(:,i) = uElems(:,i) + 0.5 * dt * (A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,1));
            else
                R = A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,i+1);
                ua(:,i) = uElems(:,i) + 0.5 * dt * (A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,i+1));
            end
        end
        for i = 1:elements
            if i == 1
                ub(:,i) = uElems(:,i) + 0.5 * dt * (A*ua(:,i) + B*ua(:,elements) + C*ua(:,i+1));
            elseif i == elements
                ub(:,i) = uElems(:,i) + 0.5 * dt * (A*ua(:,i) + B*ua(:,i-1) + C*ua(:,1));
            else
                ub(:,i) = uElems(:,i) + 0.5 * dt * (A*ua(:,i) + B*ua(:,i-1) + C*ua(:,i+1));
            end
        end
        for i = 1:elements
            if i == 1
                uc(:,i) =  uElems(:,i) + dt * (A*ub(:,i) + B*ub(:,elements) + C*ub(:,i+1));
            elseif i == elements
                uc(:,i) =  uElems(:,i) + dt * (A*ub(:,i) + B*ub(:,i-1) + C*ub(:,1));
            else
                uc(:,i) = uElems(:,i) + dt * (A*ub(:,i) + B*ub(:,i-1) + C*ub(:,i+1));
            end
        end
        for i = 1:elements
            if i == 1
                ufin(:,i) = uElems(:,i) + (1/6)*dt*( A*uElems(:,i) + B*uElems(:,elements) + C*uElems(:,i+1)   +   2*( A*ua(:,i) + B*ua(:,elements) + C*ua(:,i+1)   +   A*ub(:,i) + B*ub(:,elements) + C*ub(:,i+1) )   +    A*uc(:,i) + B*uc(:,elements) + C*uc(:,i+1)  );
            elseif i == elements
                ufin(:,i) = uElems(:,i) + (1/6)*dt*( A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,1)   +   2*( A*ua(:,i) + B*ua(:,i-1) + C*ua(:,1)   +   A*ub(:,i) + B*ub(:,i-1) + C*ub(:,1) )   +    A*uc(:,i) + B*uc(:,i-1) + C*uc(:,1)  );                
            else
                ufin(:,i) = uElems(:,i) + (1/6)*dt*( A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,i+1)   +   2*( A*ua(:,i) + B*ua(:,i-1) + C*ua(:,i+1)   +   A*ub(:,i) + B*ub(:,i-1) + C*ub(:,i+1) )   +    A*uc(:,i) + B*uc(:,i-1) + C*uc(:,i+1)  );
            end
        end
        uElems = ufin;
        uhist(:,j) = reshape(uElems, [NL*elements,1]);
    end
      
    % result -> time-marched u     
    u = reshape(uElems, [NL*elements,1]);
    Hg = zeros(NL*elements, NL*elements);
    for i=1:elements
        Hg( (i-1)*NL+1 : i*NL , (i-1)*NL+1 : i*NL) = H;
    end  
    uexact = NaN(length(x),1);
    for i=1:length(uexact)
        uexact(i) = exp( -0.5*((x(i)-0.5)/0.08)^2 );
    end
       
    % error
    err = sqrt((u-uexact)' * Hg * (u-uexact));
    
    % energy
    uHu = NaN(N,1);
    for i=1:N
        uHu(i) = uhist(:,i)' * Hg * uhist(:,i);%- uexact' * Hg * uexact;
    end    

    % t vector
    t = NaN(N,1);
    for i=1:N
        t(i) = t0 + dt*double(i);
    end
    
    % 1Transpose * H * u
    oneHu = NaN(N,1);
    for i=1:N
        oneHu(i) = ones(NL*elements,1)' * Hg * uhist(:,i);
    end   

end

        
        
        
        
        