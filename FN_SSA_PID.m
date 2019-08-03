function [t,x] = FN_SSA_PID(t,x,k)
    
    %% Perturbation
    k.(k.Pn) = k.(k.Pn)*k.P(2,sum([t>=k.P(1,:)]));
    
    %% Species:
    Z1   = x(1);
    Z2   = x(2);
    X1   = x(3);
    XC   = x(4);
    A    = x(5);
    M    = x(6);
    
    %% Functions:
    O0th = @(c)     [c];
    O1th = @(X,c)   [c*X];
    O2th = @(X,Y,c) [c*X*Y];
    OMMp = @(X,c,d) [c*X/(X+d)];
    OMMn = @(X,c,d) [c*d/(X+d)];
    
    %% Update reaction propensities:
    r = zeros(13,1);
    r(1) = O2th(Z1,Z2,k.et);
    r(2) = O0th(k.mu*k.Y);
    r(3) = O1th(Z1,k.gD);
    r(4) = O1th(XC,k.th);
    r(5) = O1th(Z2,k.gD);
    r(6) = O1th(Z1,k.bI) + OMMn(k.th*XC,k.bP*k.Y,k.mu*k.Y) + O1th(A,k.bD);
    r(7) = O1th(X1,k.g1) + O1th(X1,k.gD);
    r(8) = O1th(X1,k.bC);
    r(9) = O1th(XC,k.gC) + O1th(XC,k.gD);
    r(10)= O1th(M,k.bA);
    r(11)= OMMp(A,XC*k.gA,k.KA) + O1th(A,k.gA0) + O1th(A,k.gD);
    r(12)= O0th(k.bM*k.Y);
    r(13)= OMMp(M,A*k.gM,k.KM) + O1th(M,k.gD);
    
    %% Update time:
    % Time for the next reaction:
    tR = -log(1-rand())/sum(r);
    % Actualize global time:
    t = t + tR;
    
    %% Update system:
    syn = @(X) [X+1];
    deg = @(X) [X-1];
    rR = sum([rand()>cumsum(r)/sum(r)])+1;
    if (rR == 1)
        x(1:2) = deg(x(1:2));
    elseif (mod(rR,2)==0)
        x(rR/2) = syn(x(rR/2));
    else
        x(floor(rR/2)) = deg(x(floor(rR/2)));
    end
    