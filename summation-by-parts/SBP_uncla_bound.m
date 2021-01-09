% Element-type SBP operators with boundary nodes.
% Errorless to my knowlege.

close all; clear; clc

a = 1;
tN = 100;
t0 = 0;
x0 = 0;
xN = 1;
sigma0 = 0;
sigma1 = 1;
elements1 = 5;
elements2 = 30;
elements3 = 40;
elements4 = 50;
[u1 x1 err1 uHu1 oneHu1 t1 dx1] = SBP_solver (a, tN, t0, x0, xN, sigma1, elements1);
% [u2 x2 err2 uHu2 oneHu2 t2 dx2] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements2);
% [u3 x3 err3 uHu3 oneHu3 t3 dx3] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements3);
% [u4 x4 err4 uHu4 oneHu4 t4 dx4] = SBP_solver (a, tN, t0, x0, xN, sigma0, elements4);

% %plot
% plot(x1,u1);
% xlabel('x');
% ylabel('u');
% hold on
% xlim([0,1]);
% xzero = 0:1/50:1;
% xzero = xzero';
% uzero = NaN(length(xzero),1);
% for i=1:length(uzero)
%     uzero(i) = exp(-0.5*((xzero(i)-0.5)/0.08)^2);
% end
% % hold on
% % plot(x2,u2)
% %plot(x2,u2)
% %plot(x3,u3)
% %plot(x4,u4)
% plot(xzero, uzero);
% legend ('t=1','t=5', 'exact soln');
% title('u(x,5)');
% 
% % error - loglog plot
% deltax = [dx1 dx2 dx3 dx4]';
% errvals = [err1 err2 err3 err4]';
% loglog(deltax, errvals, 'x')
% xlabel('delta x');
% ylabel('error');
% hold on
% 
% p = polyfit(log(deltax),log(errvals),1)
% err_ = polyval(p,log(deltax));
% loglog(deltax,exp(err_));
% [m,b] = polyfit(log(deltax),log(errvals),1);

% 
% %uTranspose * H * u
% plot(t1,uHu1);
% % hold on
% % plot(t4, uHu4)
% % title('symmetric energy uTranspose*H*u vs. t');
% % legend('p = 1','p = 2');

% % onesTranspose * H * u
% %close all; clear; clc
% plot(t1,oneHu1);
% % hold on
% % plot(t4,oneHu4);
% % plot(t3,oneHu3);
% % plot(t4,oneHu4);
% % title('p=2 onesTranspose*H*u vs. t');
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
%        NL is the number of nodes per element
%        elements is the number of elements
%
%

function[u x err uHu oneHu t dx] = SBP_solver (a, tN, t0, x0, xN, sigma, elements)

    % x
    x = NaN(5,elements);
    dElems = (xN - x0)/elements;
    for i = 1:elements
        x(1,i) = x0 + dElems*(i-1);
    end
    for i = 1:elements
        if i==elements
            x(5,i) = xN;
        else
            x(5,i) = x(1,i+1);
        end
        x(2,i) = (x(5,i)-x(1,i)) / (1-(-1)) * (-sqrt(21)/7-1) + x(5,i);
        x(3,i) = (x(5,i)-x(1,i)) / (1-(-1)) * (0-1) + x(5,i);
        x(4,i) = (x(5,i)-x(1,i)) / (1-(-1)) * (sqrt(21)/7-1) + x(5,i);
    end
    xElems = x;
    x = reshape(xElems,[5*elements,1]);
      
    % u init
    uElems = NaN(5, elements);
    for i = 1:5
        for j = 1:elements
            uElems(i,j) = exp(-0.5*( (xElems(i,j) - 0.5) / 0.08 )^2);
        end
    end
    
    % Dz matrix
    Dz = [-5   49/12+7*sqrt(21)/12   -8/3   49/12-7*sqrt(21)/12   -1/2;
        
          -3/4-3*sqrt(21)/28   0   8*sqrt(21)/21   -sqrt(21)/6   3/4-3*sqrt(21)/28;
          
          3/8   -7*sqrt(21)/24   0   7*sqrt(21)/24   -3/8;
          
          -3/4+3*sqrt(21)/28   sqrt(21)/6   -8*sqrt(21)/21   0   3/4+3*sqrt(21)/28;
          
          1/2   -49/12+7*sqrt(21)/12   8/3   -49/12-7*sqrt(21)/12   5];   
    J = Dz*xElems(:,1);
    D = (1./J).*Dz;
        
    % H matrix
    v = [1/10   49/90   32/45   49/90   1/10];
    Hz = zeros(5) + diag(v,0);
    H = J.*Hz;
    
    % Q matrix
    Q = H*D;
           
    % elementary vectors 
    NL = 5;
    e1 = zeros(NL,1);
    e1(1) = 1;
    eNL = zeros(NL,1);
    eNL(NL) = 1;
    
    % case upwind    
    if sigma == 0
        
        % A matrix
        A = -a*(D + inv(H)*e1*e1');
            
        % B matrix
        B = inv(H)*a*e1*eNL';     
        
        % C matrix
        C = zeros(NL);
        
    % case symmetric  
    elseif sigma == 1
        
        % A matrix
        A = (-D*a - 0.5*inv(H)*a*e1*e1' + 0.5*inv(H)*a*eNL*eNL');
                    
        % B matrix
        B = 0.5*inv(H)*a*e1*eNL';
        
        % C matrix
        C = -0.5*inv(H)*a*eNL*e1';

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
        %dummy = 0;
        for i = 1:elements 
            if i == 1
                R = A*uElems(:,i) + B*uElems(:,elements)+ C*uElems(:,i+1);
                ua(:,i) = uElems(:,i) + 0.5 * dt * (A*uElems(:,i) + B*uElems(:,elements) + C*uElems(:,i+1));
            elseif i == elements
                R = A*uElems(:,i) + B*uElems(:,i-1) + + C*uElems(:,1);
                ua(:,i) = uElems(:,i) + 0.5 * dt * (A*uElems(:,i) + B*uElems(:,i-1) + C*uElems(:,1));
            else
                R = A*uElems(:,i) + B*uElems(:,i-1) + + C*uElems(:,i+1);
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
    
    % uTranspose * H * u
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
