%% DG using Lagrange basis functions, p >= 1

% close all; clear; clc
% 
% a = 1;
% tN = 5;
% t0 = 0;
% x0 = 0;
% xN = 1;
% Nv = 2;
% alpha0 = 0; % upwind
% alpha1 = 1; % symm
% elements1 = 5;
% elements2 = 30;
% elements3 = 40;
% elements4 = 50;
% % [u1 uexact1 x1 err1 dx1 energy1 t M1 N1 uhist1] = timemarch(elements1, Nv, xN, x0, alpha0, tN, t0, a);
% % [u2 uexact2 x2 err2 dx2 energy2 t M2 N2 uhist2] = timemarch(elements2, Nv, xN, x0, alpha0, tN, t0, a);
% % [u3 uexact3 x3 err3 dx3 energy3 t M3 N3 uhist3] = timemarch(elements3, Nv, xN, x0, alpha0, tN, t0, a);
% % [u4 uexact4 x4 err4 dx4 energy4 t M4 N4 uhist4] = timemarch(elements4, Nv, xN, x0, alpha0, tN, t0, a);
% 
% % %plot
% % figure (2)
% % plot(x1,u1);
% % xlabel('x');
% % ylabel('u');
% % hold off
% % % % plot(x1, uexact1);
% 
% % figure(1)
% % plot(x1,uexact1,'b.','LineWidth',1.5,'MarkerSize', 20);
% % xlabel('$x$','Interpreter','LaTeX','Fontsize',12);
% % ylabel('$u(t=1)$','Interpreter','LaTeX','Fontsize',12);
% % hold on
% % plot(x1, u1,'r','LineWidth',1.5);
% % hold off
% % grid on;
% % box on;
% 
% % xzero = 0:1/50:1;
% % xzero = xzero';
% % uzero = NaN(length(xzero),1);
% % for i=1:length(uzero)
% %     uzero(i) = exp(-0.5*((xzero(i)-0.5)/0.08)^2);
% % end
% % hold on
% % plot(xzero, uzero);
% % legend ('t=1','t=5', 'exact soln');
% % 
% % % error - loglog plot
% % %figure(2)
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
% % hold off
% 
% % % onesTranspose * H * u
% [oneMu t] = Conservation(elements1, Nv, xN, x0, alpha1, tN, t0, a);
% plot(t, oneMu)


% -------------------------functions begin---------------------------------

function [u uexact x err dx M N uhist dt] = test_file(elements, Nv, xN, x0, alpha, tN, t0, a)
% Purpose: time march and return main result, and plot it
% Notes: IC is fixed
%        periodic BCs
%        Nv = N+1 = order of polynomial; Nv = number of nodes per element
%        symmetric: alpha = 1 | upwind: alpha = 0

    % N
    orderN = Nv-1;

    % x
    [xElems, x] = MeshUp(orderN, elements, x0, xN);
    
    % uHat init
    uElems = NaN(Nv,elements);
    for i = 1:Nv
        for j = 1:elements
            uElems(i,j) = exp(-0.5*( (xElems(i,j) - 0.5) / 0.08 )^2);
        end
    end
    
    % lR and lL
    lR = zeros(Nv,1); lL = zeros(Nv,1);
    lR(Nv) = 1; lL(1) = 1;
    
    % matrices
    [r] = RefDomain(Nv);
    [M S] = Matrices(orderN, r);
    Dr = inv(M)*S;
    [Dr_new] = makeDr(orderN, r);
    
    % Jk
    dElems = (xN - x0)/elements;
    hk = dElems;
    Jk = hk/2;

    % Time March
    dx = dElems/(Nv-1);
    dt = 0.1*dx;
    N = floor((tN-t0)/dt);
    dt = (tN-t0)/N;
    uhist = NaN(Nv*elements, N);
    for j = 1:N
        du = NaN(Nv,elements);
        ua = NaN(Nv,elements); dua = NaN(Nv,elements);
        ub = NaN(Nv,elements); dub = NaN(Nv,elements); 
        uc = NaN(Nv,elements); duc = NaN(Nv,elements); 
        ufin = NaN(Nv,elements); 
        for i = 1:elements 
            if i == 1
                fluxL = Flux(alpha, a, uElems(Nv,elements), uElems(1,i)); fluxR = Flux(alpha, a, uElems(Nv,i), uElems(1,i+1));
                dummy = uElems(Nv,i);
            elseif i == elements
                fluxL = Flux(alpha, a, uElems(Nv,i-1), uElems(1,i)); fluxR = Flux(alpha, a, uElems(Nv,i), uElems(1,1)); 
                dummy = uElems(Nv,i);
            else
                fluxL = Flux(alpha, a, uElems(Nv,i-1), uElems(1,i)); fluxR = Flux(alpha, a, uElems(Nv,i), uElems(1,i+1));
                dummy = uElems(Nv,i);
            end
            %dummy;
            du(:,i) =  inv(Jk*M)*a*(-S + lR*lR' - lL*lL')*uElems(:,i) + inv(Jk*M)*(lL*fluxL - lR*fluxR); 
            ua(:,i) = uElems(:,i) + 0.5 * dt * ( du(:,i) );
        end
        for i = 1:elements
            if i == 1
                fluxL = Flux(alpha, a, ua(Nv,elements), ua(1,i)); fluxR = Flux(alpha, a, ua(Nv,i), ua(1,i+1));
            elseif i == elements
                fluxL = Flux(alpha, a, ua(Nv,i-1), ua(1,i)); fluxR = Flux(alpha, a, ua(Nv,i), ua(1,1));
            else
                fluxL = Flux(alpha, a, ua(Nv,i-1), ua(1,i)); fluxR = Flux(alpha, a, ua(Nv,i), ua(1,i+1));
            end
            dua(:,i) =  inv(Jk*M)*a*(-S + lR*lR' - lL*lL')*ua(:,i) + inv(Jk*M)*(lL*fluxL - lR*fluxR); 
            ub(:,i) = uElems(:,i) + 0.5 * dt * ( dua(:,i) );
        end
        for i = 1:elements
            if i == 1
                fluxL = Flux(alpha, a, ub(Nv,elements), ub(1,i)); fluxR = Flux(alpha, a, ub(Nv,i), ub(1,i+1));
            elseif i == elements
                fluxL = Flux(alpha, a, ub(Nv,i-1), ub(1,i)); fluxR = Flux(alpha, a, ub(Nv,i), ub(1,1));
            else
                fluxL = Flux(alpha, a, ub(Nv,i-1), ub(1,i)); fluxR = Flux(alpha, a, ub(Nv,i), ub(1,i+1));
            end
            dub(:,i) = inv(Jk*M)*a*(-S + lR*lR' - lL*lL')*ub(:,i) + inv(Jk*M)*(lL*fluxL - lR*fluxR);
            uc(:,i) = uElems(:,i) + dt * ( dub(:,i) );
        end
        for i = 1:elements
            if i == 1
                fluxL = Flux(alpha, a, uc(Nv,elements), uc(1,i)); fluxR = Flux(alpha, a, uc(Nv,i), uc(1,i+1));
            elseif i == elements
                fluxL = Flux(alpha, a, uc(Nv,i-1), uc(1,i)); fluxR = Flux(alpha, a, uc(Nv,i), uc(1,1));
            else
                fluxL = Flux(alpha, a, uc(Nv,i-1), uc(1,i)); fluxR = Flux(alpha, a, uc(Nv,i), uc(1,i+1));
            end
            duc(:,i) = inv(Jk*M)*a*(-S + lR*lR' - lL*lL')*uc(:,i) + inv(Jk*M)*(lL*fluxL - lR*fluxR);
            ufin(:,i) = uElems(:,i) + (1/6) * dt * ( du(:,i) + 2*( dua(:,i) + dub(:,i) ) + duc(:,i) );
        end
        uElems = ufin;
        uhist(:,j) = reshape(uElems, [Nv*elements,1]);
    end
    
    % result -> time-marched u   
    u = reshape(uElems, [Nv*elements,1]);
    uexact = NaN(length(x),1);
    for i=1:length(uexact)
        uexact(i) = exp( -0.5*((x(i)-0.5)/0.08)^2 );
    end
        
    % error
    difference = u - uexact;
    err = 0;
    for i = 1:length(u)
        err = err + (difference(i))^2;
    end
    err = sqrt(hk*err);

end


function [l, dl] = lagrangeP(j, r, xi)
% returns val of jth lagrangeP and its derivative at point xi
% lagrange basis for ref element [-1,1]
    
    % construct syms lagrangeP
    l = 1;
    Nv = length(r);
    syms x;
    for m = 1:Nv
        if m==j
            l = l;
        else
            l = l*(x-r(m))/(r(j)-r(m));
        end
    end
    
    % construct syms differentiation of lagrangeP
    y = diff(l,x,1);
    
    % evaluate lagrangeP at point xi
    l = eval(subs(l, xi));
    
    % evaluate diff_lagrangeP at point xi
    dl = eval(subs(y, xi));
    
end

function [l] = lagrangeVector(r, xi)
% Purpose: constructs the value of lagrange polynomial basis at a point xi
% which is in [-1,1].
% Returns as a vector, l.
    l = NaN(length(r),1);
    for i = 1:length(l)
        l(i)= lagrangeP(i, r, xi);
    end
end

function [M S] = Matrices(N, r)
% Purpose: construct matrices on ref elem
    
    Nv = N+1;
    M = NaN(Nv,Nv);
    Dr = NaN(Nv, Nv);
    [r w P] = lglnodes(N);
    [r2,w2] = lglnodes(N+1);
    % M matrix
    for i = 1:Nv
        for j = 1:Nv
            M(i,j) = 0;
            S(i,j) = 0;
            for k = 1:length(w2)
                 [li, dli] = lagrangeP(i, r, r2(k));
                 [lj, dlj] = lagrangeP(j, r, r2(k));
                 M(i,j) = M(i,j) + w2(k)*li*lj;
                 S(i,j) = S(i,j) + w2(k)*li*dlj;
            end
        end
    end
end

function [Dr] = makeDr(N, r)
    Nv = N+1;
    Dr = NaN(Nv);
    [r w P] = lglnodes(N);
    for i = 1:Nv
        for j = 1:Nv
            [lj, dlj] = lagrangeP(j, r, r(i));
            Dr(i,j) = dlj;
        end
    end
end

function [x w P] = lglnodes(N)
% Purpose: computes zeros of (1-x^2)*P'_N(x)
% Creds to Greg von Winckel 

 % Truncation + 1
 N1=N+1;

 % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
 x = cos(pi*(0:N)/N)';

 % The Legendre Vandermonde Matrix
 P = zeros(N1,N1);

 % Compute P_(N) using the recursion relation
 % Compute its first and second derivatives and 
 % update x using the Newton-Raphson method.
 xold=2;
 while max(abs(x-xold))>eps
     xold=x;
     P(:,1)=1;    P(:,2)=x;
     for k=2: N
         P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
     end
     x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
 end
 w=2./(N*N1*P(:,N1).^2);
 x = flipud(x);
end

function [r] = RefDomain(Nv)
% Purpose: compute the reference domain vector given the number of nodes
% per element, Nv
    [r,w,P] = lglnodes(Nv-1);
end

function [xElems x] = MeshUp(N, elements, x0, xN) 
% Purpose: returns the mesh
    [xi] = RefDomain(N+1); Nv = N+1; % I use xi = r
    x = NaN(Nv,elements);
    xmax = NaN(elements,1);
    xmin = NaN(elements,1);
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
        x(:,i) = (xmax(i)-xmin(i)) / (1-(-1)) * (xi-1) +  xmax(i);
    end
    xElems = x;
    x = reshape(xElems,[Nv*elements,1]);
end

function [flux] = Flux(alpha, a, uL, uR)
% Purpose: returns flux as a scalar
% symmetric: alpha = 1
% upwind:    alpha = 0

    flux = a*(uL+uR)/2 + abs(a)*((1-alpha)/2)*(uL-uR);
    
end

function [energy t] = Energy(elements, Nv, xN, x0, alpha, tN, t0, a)
% Purpose: returns energy vs time and plots it
    
    % call time march function
    [u uexact x err dx M N uhist dt] = timemarch(elements, Nv, xN, x0, alpha, tN, t0, a);
    
    % Jk
    dElems = (xN - x0)/elements;
    hk = dElems;
    Jk = hk/2;
    
    % energy
    Mg = zeros(Nv*elements, Nv*elements); % global mass matrix
    Mk = Jk*M;
    for i=1:elements
        Mg( (i-1)*Nv+1 : i*Nv , (i-1)*Nv+1 : i*Nv) = Mk;
    end  
    energy = NaN(N,1);
    for i=1:N
        energy(i) = uhist(:,i)' * Mg * uhist(:,i);
    end    

    % t vector
    t = NaN(N,1);
    for i=1:N
        t(i) = t0 + dt*double(i);
    end
    
end

function [oneMu t] = Conservation(elements, Nv, xN, x0, alpha, tN, t0, a)
    
    % Call time march function
    [u uexact x err dx M N uhist dt] = timemarch(elements, Nv, xN, x0, alpha, tN, t0, a);
    
    % Jk
    dElems = (xN - x0)/elements;
    hk = dElems;
    Jk = hk/2;
    
    % Conservation
    Mg = zeros(Nv*elements, Nv*elements); % global mass matrix
    Mk = Jk*M;
    
    oneMu = NaN(N,1);
    for i=1:N
        oneMu(i) = ones(Nv*elements,1)' * Mg * uhist(:,i);
    end
    
    % t vector
    t = NaN(N,1);
    for i=1:N
        t(i) = t0 + dt*double(i);
    end
end
    

    





