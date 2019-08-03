%% ODE linearized model of the PID controller designed in Chevalier, Gómez-Schiavon et al. (2019; Cell Systems)
function dYdt = FN_ODE_PID_Lin(t,y)
    load Par_ODE.mat
	
	Ys  = k.Y;
    
    %% Perturbation
    k.(k.Pn) = k.(k.Pn)*k.P(2,sum([t>=k.P(1,:)]));
	
	Y   = k.Y - Ys;

    %% Species:
    z1   = y(1);
    z2   = y(2);
    x1   = y(3);
    xC   = y(4);
    a    = y(5);
    m    = y(6);

    %% Steady states:
    Z1s  = k.SS(1);
    Z2s  = k.SS(2);
    X1s  = k.SS(3);
    XCs  = k.SS(4);
    As   = k.SS(5);
    Ms   = k.SS(6);
    
    %% Events/functions:
    O0th = @(c)     [c];
    O1th = @(X,c)   [c*X];
    O2th = @(X,Y,c) [c*X*Y];
    OMMp = @(X,c,d) [c*X/(X+d)];
    OMMn = @(X,c,d) [c*d/(X+d)];
    
    OMLn = @(X,c,d) [c*d/((X+d)^2)];
    
    %% ODEs
    F0 = [O0th(k.mu*Ys)  - O2th(Z1s,Z2s,k.et)    - O1th(Z1s,k.gD);
          O1th(XCs,k.th) - O2th(Z1s,Z2s,k.et)    - O1th(Z2s,k.gD);
          O1th(Z1s,k.bI)  ...
            + OMMn(k.th*XCs,k.bP*Ys,k.mu*Ys) + O1th(As,k.bD) ...
                             - O1th(X1s,k.g1)    - O1th(X1s,k.gD);
          O1th(X1s,k.bC)  - O1th(XCs,k.gC)       - O1th(XCs,k.gD);
          O1th(Ms ,k.bA)  - OMMp(As,XCs*k.gA,k.KA) ...
                             - O1th(As,k.gA0)    - O1th(As,k.gD);
          O0th(k.bM*Ys)  - OMMp(Ms,As*k.gM,k.KM) - O1th(Ms,k.gD)];
			
    dYdt = F0 + ...
	       [O0th(k.mu*Y)   - O1th(z1,k.et*Z2s) - O1th(z2,k.et*Z1s) ...
												 - O1th(z1,k.gD);
            O1th(xC,k.th)  - O1th(z1,k.et*Z2s) - O1th(z2,k.et*Z1s) ...
												 - O1th(z2,k.gD);
            O1th(z1,k.bI)  ...
              + (k.bP*Y*(OMMn(k.th*XCs,2      ,k.mu*Ys) ...
			           - OMLn(k.th*XCs,k.mu*Ys,k.mu*Ys)))...
			  - (k.th*xC*OMLn(k.th*XCs,k.bP*Ys,k.mu*Ys))...
			  + O1th(a,k.bD) - O1th(x1,k.g1)    - O1th(x1,k.gD);
            O1th(x1,k.bC)  - O1th(xC,k.gC)       - O1th(xC,k.gD);
            O1th(m ,k.bA) ...
			  - OMMp(As,xC*k.gA,k.KA) - OMLn(As,a*XCs*k.gA,k.KA) ...
                              - O1th(a,k.gA0)    - O1th(a ,k.gD);
            O0th(k.bM*Y) ...
			  - OMMp(Ms,a*k.gM,k.KM) - OMLn(Ms,m*As*k.gM,k.KM) ...
												 - O1th(m ,k.gD)];
    