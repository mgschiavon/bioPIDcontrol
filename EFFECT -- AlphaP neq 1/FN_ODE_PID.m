function dYdt = FN_ODE_PID(t,y)
    load Par_ODE.mat
    
    if(y(2)>4000)
        myOut = 'Z_2 too large'
        y = y * NaN;
    end
    
    %% Perturbation
    k.(k.Pn) = k.(k.Pn)*k.P(2,sum([t>=k.P(1,:)]));

    %% Species:
    Z1   = y(1);
    Z2   = y(2);
    X1   = y(3);
    XC   = y(4);
    A    = y(5);
    M    = y(6);
    
    %% Events/functions:
    O0th = @(c)     [c];
    O1th = @(X,c)   [c*X];
    O2th = @(X,Y,c) [c*X*Y];
    OMMp = @(X,c,d) [c*X/(X+d)];
    OMMn = @(X,c,d) [c*d/(X+d)];
    
    %% ODEs
    dYdt = [O0th(k.mu*k.Y) - O2th(Z1,Z2,k.et)    - O1th(Z1,k.gD);
            O1th(XC,k.th)  - O2th(Z1,Z2,k.et)    - O1th(Z2,k.gD);
            O1th(Z1,k.bI)  ...
            + OMMn(k.th*XC,k.bP*k.Y,k.aP*k.mu*k.Y) + O1th(A,k.bD) ...
                              - O1th(X1,k.g1)    - O1th(X1,k.gD);
            O1th(X1,k.bC)  - O1th(XC,k.gC)       - O1th(XC,k.gD);
            O1th(M ,k.bA)  - OMMp(A,XC*k.gA,k.KA) ...
                              - O1th(A,k.gA0)    - O1th(A,k.gD);
            O0th(k.bM*k.Y) - OMMp(M,A*k.gM,k.KM) - O1th(M,k.gD)];
    