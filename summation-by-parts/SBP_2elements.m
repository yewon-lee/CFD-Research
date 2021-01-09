% This is a fully functional (errorless) version of the SBP program with
% two elements only.

close all; clear; %clc

a = 1;
p1 = 1;
p2 = 2;
NL1 = 50;
NL2 = 100;
NL3 = 200;
NL4 = 500;
tN = 1;
t0 = 0;
x0 = 0;
xN = 1;
sigma = 0;
[u1 x1 err1 uHu1 oneHu1 t1 ] = SBP_solver (a, p1, NL1, tN, t0, x0, xN, sigma);
% [u2 x2 err2 uHu2 oneHu2 t2 ] = SBP_solver (a, p2, NL2, tN, t0, x0, xN, sigma);
% [u3 x3 err3 uHu3 oneHu3 t3 ] = SBP_solver (a, p2, NL3, tN, t0, x0, xN, sigma);
% [u4 x4 err4 uHu4 oneHu4 t4 ] = SBP_solver (a, p2, NL4, tN, t0, x0, xN, sigma);

% plot
% plot(x1,u1);
% xlabel('x');
% ylabel('u');
% hold on
plot(x1,u1);
xzero = 0:1/50:1;
xzero = xzero';
uzero = NaN(length(xzero),1);
for i=1:length(uzero)
    uzero(i) = exp(-0.5*((xzero(i)-0.5)/0.08)^2);
end
hold on
%plot(x2,u2)
plot(xzero, uzero);
legend ('p=1', 'exact soln');
title('u(x,5)');
% 
% % error - loglog plot
% deltax = [1/NL1 1/NL2 1/NL3 1/NL4]';
% errvals = [err1 err2 err3 err4]';
% loglog(deltax, errvals, 'x')
% xlabel('delta x');
% ylabel('error');
% hold on
% 
% p = polyfit(log(deltax),log(errvals),1)
% err_ = polyval(p,log(deltax));
% loglog(deltax,exp(err_));
% % [m,b] = polyfit(log(deltax),log(errvals),1);
% % disp(m);

% % uTranspose * H * u
% %close all; clear; clc
% plot(t1,uHu1);
% ylim([0.15 0.25]);
% hold on
% plot(t2,uHu2);
% title('uTranspose*H*u vs. t');
% legend('p = 1','p = 2');
% A = [t1,uHu1];

% % onesTranspose * H * u
% %close all; clear; clc
% plot(t1,oneHu1);
% hold on
% plot(t2,oneHu2);
% title('onesTranspose*H*u vs. t');
% legend('p = 1','p = 2');

% ----------------------function---------------------

% Returns u(x) at t=tN
%
% Returns rate of convergence with grid refinement of the H-norm of the
% solution error
%
% Returns uTranspose * H * u as a function of time
%
% Returns 1Transpose * H * u as a function of time
%
% Note: sigma = 1 -> symmetric
%       sigma = 0 -> upwind
%
%

function[u x err uHu oneHu t] = SBP_solver (a, p, NL, tN, t0, x0, xN, sigma)

    % u init
    ufun = @(x) exp(-0.5*((x-0.5)/0.08)^2);
    NR = NL;
    uL = NaN(NL,1);
    xM = (xN-x0)/2;
    dxL = (xM-x0)/(NL-1);
    dxR = (xN-xM)/(NR-1);
    uR = NaN(NR,1);
    for i=1:NL
        uL(i) = ufun(x0 + dxL*(i-1));
    end
    for i=1:NR
        uR(i) = ufun(xM + dxR*(i-1));
    end
    
    % x
    xL = linspace(x0,xM,NL)';
    xR = linspace(xM,xN,NL)';
    x = [xL;xR];
    
    % case of p=1
    if p == 1
        
    % left
    
        h = dxL;
        
        % DL matrix
        DL = (1/h) * (zeros(NL,NL) + diag((-1/2)*ones(NL-1,1),-1) + diag((1/2)*ones(NL-1,1),1) );
        DL(1,1) = -1/h;
        DL(1,2) = 1/h;
        DL(NL,NL-1) = -1/h;
        DL(NL,NL) = 1/h;        
        
        % HL matrix
        v = ones(NL,1);
        v(1) = 0.5;
        v(NL) = 0.5;
        HL = h * (zeros(NL) + diag(v,0));
        
       
    % right
        
        h = dxR;
        
        % DR matrix
        DR = (1/h) * (zeros(NR) + diag((-1/2)*ones(NR-1,1),-1) + diag((1/2)*ones(NR-1,1),1) );
        DR(1,1) = -1/h;
        DR(1,2) = 1/h;
        DR(NR,NR-1) = -1/h;
        DR(NR,NR) = 1/h; 
        
        % HR matrix
        v = ones(NR,1);
        v(1) = 0.5;
        v(NR) = 0.5;
        HR = h * (zeros(NR) + diag(v,0));
        
    % case of p=2
   
    elseif p == 2
        
    % left
    
        h = dxL;
        
        % HL matrix
        v = ones(NL,1);
        v(1,1) = 17/48;
        v(2,1) = 59/48;
        v(3,1) = 43/48;
        v(4,1) = 49/48;
        for i = 0:3
            v(NL-i) = v(i+1,1);
        end
        HL = zeros(NL) + diag(v,0);
        HL = h*HL;
        
        % QL matrix
        QL = zeros(NL) + diag((-8/12)*ones(NL-1,1),-1) + diag((8/12)*ones(NL-1,1),1) + diag((1/12)*ones(NL-2,1),-2) + diag((-1/12)*ones(NL-2,1),2);
        
        QL(1,1) = -1/2;
        QL (1,2) = 59/96;
        QL(1,3) = -1/12;
        QL(1,4) = -1/32;
        QL(2,1) = -59/96;
        QL(2,2) = 0;
        QL(2,3) = 59/96;
        QL(2,4) = 0;
        QL(3,1) = 1/12;
        QL(3,2) = -59/96;
        QL(3,3) = 0;
        QL(3,4) = 59/96;
        QL(4,1) = 1/32;
        QL(4,2) = 0;
        QL(4,3) = -59/96;
        QL(4,4) = 0;
        
        QL(NL-3,NL-3) = 0;
        QL(NL-3,NL-2) = 59/96;
        QL(NL-3,NL-1) = 0;
        QL(NL-3,NL) = -1/32;
        QL(NL-2,NL-3) = -59/96;
        QL(NL-2,NL-2) = 0;
        QL(NL-2,NL-1) = 59/96;
        QL(NL-2,NL) = -1/12; 
        QL(NL-1,NL-3) = 0;
        QL(NL-1,NL-2) = -59/96;
        QL(NL-1,NL-1) = 0;
        QL(NL-1,NL) = 59/96;
        QL(NL,NL-3) = 1/32;
        QL(NL,NL-2) = 1/12;
        QL(NL,NL-1) = -59/96;
        QL(NL,NL) = 1/2;
        
        % DL matrix
        DL = inv(HL)*QL;
      
    % right
        
        h = dxR;
        
        % HR matrix
        v = ones(NR,1);
        v(1,1) = 17/48;
        v(2,1) = 59/48;
        v(3,1) = 43/48;
        v(4,1) = 49/48;
        for i = 0:3
            v(NR-i) = v(i+1,1);
        end
        HR = zeros(NR) + diag(v,0);
        HR = h*HR;
        
        % QR matrix
        QR = QL;
        
        % DR matrix
        DR = inv(HR)*QR;        
        
    end
                
           
    % case upwind
    
    if sigma == 0
        
        % AL matrix
        e1L = zeros(NL,1);
        e1L(1) = 1;
        eNR = zeros(NR,1);
        eNR(NR) = 1;
        AL = -a*(DL + inv(HL)*e1L*e1L');
            
        % BL matrix
        BL = inv(HL)*a*e1L*eNR';
    
        % AR matrix
        e1R = zeros(NR,1);
        e1R(1) = 1;
        eNL = zeros(NL,1);
        eNL(NL) = 1;
        AR = -a*(DR + inv(HR)*e1R*e1R');
            
        % BR matrix
        BR = inv(HR)*a*e1R*eNL';
        
        
    % case symmetric  
    elseif sigma == 1
        
        % AL matrix
        e1L = zeros(NL,1);
        e1L(1) = 1;
        eNR = zeros(NR,1);
        eNR(NR) = 1;
        AL = (-DL*a - 0.5*inv(HL)*a*e1L*e1L');
                    
        % BL matrix
        BL = 0.5*inv(HL)*a*e1L*eNR';
    
        % AR matrix
        e1R = zeros(NR,1);
        e1R(1) = 1;
        eNL = zeros(NL,1);
        eNL(NL) = 1;
        AR = -a*(DR + 0.5*inv(HR)*e1R*e1R');
            
        % BR matrix
        BR = 0.5*inv(HR)*a*e1R*eNL';
    end
        
    % time marching
    dt = 0.1*dxL;
    N = (tN-t0)/dt;
    N = floor((tN-t0)/dt);
    dt = (tN-t0)/N;
    uMatrix = NaN(length(uL) + length(uR), N);
    for i = 1:N
        uaL = uL + 0.5 * dt * (AL*uL + BL*uR);
        uaR = uR + 0.5 * dt * (AR*uR + BR*uL);
        ubL = uL + 0.5 * dt * (AL*uaL + BL*uaR);
        ubR = uR + 0.5 * dt * (AR*uaR + BR*uaL);
        ucL = uL + dt * (AL*ubL + BL*ubR);
        ucR = uR + dt * (AR*ubR + BR*ubL);
        uLfin = uL + (1/6)*dt*( AL*uL + BL*uR   +   2*( AL*uaL + BL*uaR   +   AL*ubL + BL*ubR)   +    AL*ucL + BL*ucR  );
        uRfin = uR + (1/6)*dt*( AR*uR + BR*uL   +   2*( AR*uaR + BR*uaL   +   AR*ubR + BR*ubL)   +    AR*ucR + BR*ucL  ); 
        uL = uLfin;
        uR = uRfin;
        uMatrix(:,i) = [uL;uR];  
    end
    
    % result -> time-marched u     
    u = [uL;uR];
    H = blkdiag(HL,HR);
    uexact = NaN(length(u),1);
    for i=1:length(uexact)
        %uexact(i) = exp( -0.5*((x(i)-0.5)/0.08)^2 );
        uexact(i) = ufun(x(i));
    end
       
    % error
    err = sqrt((u-uexact)' * H * (u-uexact));
    
    % uTranspose * H * u
    uHu = NaN(N,1);
    for i=1:N
        uHu(i) = uMatrix(:,i)' * H * uMatrix(:,i);
    end
    
    % t vector
    t = NaN(N,1);
    for i=1:N
        t(i) = t0 + dt*double(i-1);
    end
    
    % 1Transpose * H * u
    oneHu = NaN(N,1);
    for i=1:N
        uHu(i) = ones(length(u),1)' * H * uMatrix(:,i);
    end
end

        
        