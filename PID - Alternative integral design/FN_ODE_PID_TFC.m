%% ODE system of the PID controller designed in Chevalier, Gómez-Schiavon et al. (2019; Cell Systems)
%  using a new integral controller design instead of the antithetic motif
function dYdt = FN_ODE_PID_TFC(t,y)
    load Par_ODE.mat
    
    %% Perturbation
    k.(k.Pn) = k.(k.Pn)*k.P(2,sum([t>=k.P(1,:)]));

    %% Species:
    X1   = y(1);
    A    = y(2);
    M    = y(3);
    Z    = y(4);
    XC   = y(5);
    
    %% Events/functions:
    O0th  = @(c)     [c];
    O1th  = @(X,c)   [c*X];
    OMMp  = @(X,c,d) [c*X/(X+d)];
    OMMn  = @(X,c,d) [c*d/(X+d)];
    
    %% ODEs
    dYdt = [O1th(Z,k.bI) + OMMn(k.th*XC,k.bP*k.Y,k.mu*k.Y)...
                         + O1th(A,k.bD)...
                            - O1th(X1,k.g1) - O1th(X1,k.g);    % X1
            O1th(M,k.bA) - OMMp(A,XC*k.gA,k.KA)...
                            - O1th(A,k.gA0) - O1th(A ,k.g);    % A
            O0th(k.bM*k.Y) - OMMp(M,A*k.gM,k.KM)...
                                            - O1th(M ,k.g);    % M
            OMMn(k.th*XC,k.bZ*k.Y,k.mu*k.Y)...
                           - k.gZ*k.aZ*k.Y*Z/(Z+k.KZ) - O1th(Z,k.g); % Z
            O1th(X1,k.bC)   - O1th(XC,k.gC) - O1th(XC,k.g)];   % XC
    
    if ((Z+dYdt(4)) < 0)
        Z = 0;
        dYdt(4) = 0;
    end;
    