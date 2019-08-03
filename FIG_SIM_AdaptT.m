%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%% Quantifying performance %%%%%%%%
%%%%%%%%%%  Adaptation Time  %%%%%%%%%%%
% Mariana Gómez-Schiavon
% December, 2017

clear;

%% Quantification parameters:
T   = 5000;     % Maximum simulation time
eps = 0.05;     % Adaptation threshold

k.Y  = 300;
k.bI = 0.06;

sweaP1 = 'bP';  % Parameter to vary in y-axis
sweaI1 = 25;    % Number of bins
sweaM1 = 14/30; % Maximum value (e.g. such that kP = 0.035, then k.bP = 0.035*4*k.mu/k.th)
sweaP2 = 'bD';	% Parameter to vary in x-axis
sweaI2 = 25;    % Number of bins
sweaM2 = 0.9;	% Maximum value (e.g. such that kD = 0.6, then k.bD = 0.6*(k.bA*k.gM/k.gA)

%% Kinetic parameters:
k.gD = 0.01;    % Dilution
k.mu = 1;
k.et = 0.01;
k.th = 0.3;
k.g1 = 0.1;
k.bC = 0.1;
k.gC = 0.1;
k.bA = 1.5;
k.gA = 1.5;
k.KA = 1;
k.gA0= 0.1;
k.bM = 5/12;
k.gM = 1.5;
k.KM = 1;

% Perturbation
k.P  = [0;      % Time points
        1.5];   % Multiplicative perturbation
k.Pn = 'bC';    % Variable to be perturbed

%% Sweaping:
sweaS1 = sweaM1/sweaI1;
sweaS2 = sweaM2/sweaI1;
AT = zeros(sweaI1,sweaI2,2) + Inf;
mE = zeros(sweaI1,sweaI2,2) - Inf;
mY = zeros(sweaI1,sweaI2,6) - Inf;
O = k;
for i = 1:sweaI1
    for j = 1:sweaI2
        [i,j]
        k.(sweaP1) = i*sweaS1;
        k.(sweaP2) = j*sweaS2;
        % Steady states ignoring dilution
        XC = k.mu*k.Y/k.th;
        X1 = XC*k.gC/k.bC;
        A  = k.bM*k.Y/k.gM;
        M  = ((k.gA*XC)+(k.gA0*A))/k.bA;
        Z1 = ((k.g1*X1)-(k.bP*k.Y*k.mu*k.Y/((k.th*XC)...
            +(k.mu*k.Y)))-(k.bD*A))/k.bI;
        Z2 = k.th*XC/(k.et*Z1);
        SS = [Z1 Z2 X1 XC A M];
        clear Z1 Z2 X1 XC A M
        
        % ODE:
        save('Par_ODE.mat','k');
        [t,y] = ode45(@FN_ODE_PID,[0:T],SS);
        XC  = y(:,4);
        mY(i,j,:) = mean(y);
        mE(i,j,:) = [min(XC-mY(i,j,4)) max(XC-mY(i,j,4))];

        % Adaptation performance:
        eX = abs(XC-mY(i,j,4))/mY(i,j,4);
        AT(i,j,1) = sum([cumsum([eX([T:-1:0]+1)<=eps]) == [1:(T+1)]']);
        eX = abs(XC-mY(i,j,4))/max(abs(mE(i,j,:)));
        AT(i,j,2) = sum([cumsum([eX([T:-1:0]+1)<=eps]) == [1:(T+1)]']);
        clear t y SS XC eX
    end
    save(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))
end
k = O;
clear i j O ans
AT = T + 1 - AT;
save(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))

%% Heat-map figure
fig = figure();
fig.Units = 'inches';
fig.Position = [5 2 3 2];
fig.PaperPosition = fig.Position;
C = colormap('jet');
C(1,:)  = [1 1 1];
C(64,:) = [0 0 0];
fig.Colormap = C;
    at(:,:) = AT(:,:,2);
    mX(:,:) = mY(:,:,4);
    at([mX<(0.9*k.mu*k.Y/k.th)]) = NaN;
    at([mX>(1.1*k.mu*k.Y/k.th)]) = NaN;
    
    imagesc([1:sweaI2]*sweaS2*k.gA/(k.bA*k.gM),...
        [1:sweaI1]*sweaS1*k.th/(4*k.mu),at,[0 300])
        xlabel('k_D \theta','FontSize',12)
        ylabel('k_P \theta','FontSize',12)
        set(gca,'YDir','normal','FontSize',12,...
            'XTick',[0.2:0.2:0.6])
        colorbar;
clear a b c

%% END