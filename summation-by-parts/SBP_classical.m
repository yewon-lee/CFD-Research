%% SBP using classical operators, p=1 and p=2 (Q1.5)

% This is the SBP program for Question 1.5: period BC's, p=1 and 2, and RK4
% time marching. 
% 
% No errors present to my knowledge. 
%
% Latest version.
%

close all; clear; clc

a = 1;
p1 = 1;
p2 = 2;
NL1 = 5;
NL2 = 5;
NL3 = 5;
NL4 = 5;
tN = 5;
tN2 = 1;
t0 = 0;
x0 = 0;
xN = 1;
sigma0 = 0; % upwind
sigma1 = 1; % symm
elements1 = 50;
elements2 = 60;
elements3 = 70;
elements4 = 80;
% elements1 = 2^6; 
% elements2 = 2^7; elements3 = 2^8; elements4 = 2^9;
[u1 x1 err1 uHu1 oneHu1 t1 dx1] = SBP_solver (a, p2, NL1, tN, t0, x0, xN, sigma0, elements1);
% [u2 x2 err2 uHu2 oneHu2 t2 dx2] = SBP_solver (a, p2, NL1, tN, t0, x0, xN, sigma0, elements2);
% [u3 x3 err3 uHu3 oneHu3 t3 dx3] = SBP_solver (a, p2, NL1, tN, t0, x0, xN, sigma0, elements3);
% [u4 x4 err4 uHu4 oneHu4 t4 dx4] = SBP_solver (a, p2, NL1, tN, t0, x0, xN, sigma0, elements4);
% 
%plot
plot(x1,u1);
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
% hold on
% plot(x2,u2)
% %plot(x2,u2)
% %plot(x3,u3)
% %plot(x4,u4)
% plot(xzero, uzero);
% legend ('t=1','t=5', 'exact soln');
% title('u(x,5)');

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
% %disp(m);

% %uTranspose * H * u
% plot(t1,uHu1);
% 
% % % onesTranspose * H * u
% % %close all; clear; clc
% % plot(t3,oneHu3);
% % hold on
% % plot(t4,oneHu4);
% % % plot(t3,oneHu3);
% % % plot(t4,oneHu4);
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

function[u x err uHu oneHu t dx] = SBP_solver (a, p, NL, tN, t0, x0, xN, sigma, elements)

    % x
    dx = (xN-x0)/(elements)/(NL-1);
    x = NaN(NL*elements,1);
    x(1) = x0;
    accum = 0;
    for i = 2:NL*elements
        if accum == NL-1
            x(i) = x(i-1);
            accum = 0;
        else 
            x(i) = x(i-1) + dx;
            accum = accum + 1;
        end
    end
    xElems = reshape(x,[NL,elements]);
    
    % u init
    uElems = NaN(NL, elements);
    for i = 1:NL
        for j = 1:elements
            uElems(i,j) = exp(-0.5*( (xElems(i,j) - 0.5) / 0.08 )^2);
        end
    end
    
    % case of p=1
    if p == 1
       
        h = dx;
        
        % D matrix
        D = (1/h) * (zeros(NL,NL) + diag((-1/2)*ones(NL-1,1),-1) + diag((1/2)*ones(NL-1,1),1) );
        D(1,1) = -1/h;
        D(1,2) = 1/h;
        D(NL,NL-1) = -1/h;
        D(NL,NL) = 1/h;        
        
        % H matrix
        v = ones(NL,1);
        v(1) = 0.5;
        v(NL) = 0.5;
        H = h * (zeros(NL) + diag(v,0));
        
       
    % case of p=2
   
    elseif p == 2
  
        h = dx;
        
        % H matrix
        v = ones(NL,1);
        v(1,1) = 17/48;
        v(2,1) = 59/48;
        v(3,1) = 43/48;
        v(4,1) = 49/48;
        for i = 0:3
            v(NL-i) = v(i+1,1);
        end
        H = zeros(NL) + diag(v,0);
        H = h*H;
        
        % Q matrix
        Q = zeros(NL) + diag((-8/12)*ones(NL-1,1),-1) + diag((8/12)*ones(NL-1,1),1) + diag((1/12)*ones(NL-2,1),-2) + diag((-1/12)*ones(NL-2,1),2);
        
        Q(1,1) = -1/2;
        Q(1,2) = 59/96;
        Q(1,3) = -1/12;
        Q(1,4) = -1/32;
        Q(2,1) = -59/96;
        Q(2,2) = 0;
        Q(2,3) = 59/96;
        Q(2,4) = 0;
        Q(3,1) = 1/12;
        Q(3,2) = -59/96;
        Q(3,3) = 0;
        Q(3,4) = 59/96;
        Q(4,1) = 1/32;
        Q(4,2) = 0;
        Q(4,3) = -59/96;
        Q(4,4) = 0;
        
        Q(NL-3,NL-3) = 0;
        Q(NL-3,NL-2) = 59/96;
        Q(NL-3,NL-1) = 0;
        Q(NL-3,NL) = -1/32;
        Q(NL-2,NL-3) = -59/96;
        Q(NL-2,NL-2) = 0;
        Q(NL-2,NL-1) = 59/96;
        Q(NL-2,NL) = -1/12; 
        Q(NL-1,NL-3) = 0;
        Q(NL-1,NL-2) = -59/96;
        Q(NL-1,NL-1) = 0;
        Q(NL-1,NL) = 59/96;
        Q(NL,NL-3) = 1/32;
        Q(NL,NL-2) = 1/12;
        Q(NL,NL-1) = -59/96;
        Q(NL,NL) = 1/2;
        
        % D matrix
        D = inv(H)*Q;
      
    end
           
    % elementary vectors       
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
    dt = 0.1*dx;
    N = floor((tN-t0)/dt);
    dt = (tN-t0)/N;
    uhist = NaN(NL*elements, N);
    ua = NaN(NL,elements); 
    ub = NaN(NL,elements); 
    uc = NaN(NL,elements); 
    ufin = NaN(NL,elements); 
    for j = 1:N
        dummy = 0;
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
    
    % energy
    uHu = NaN(N,1);
    for i=1:N
        uHu(i) = uhist(:,i)' * Hg * uhist(:,i);
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
