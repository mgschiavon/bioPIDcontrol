%%       FIGURE : PID-control         %%
%%%%%%%%% Complex plant design %%%%%%%%%
%%%%%%%% Adaptation time sweep %%%%%%%%%
% Mariana Gómez-Schiavon
% June, 2018

clear;

%% Quantification parameters:
T   = 5000;     % Maximum simulation time
eps = 0.05;     % Adaptation threshold

k.Y  = 600;
k.bI = 0.02;

% kP varying:
sIkP = [0:(0.14/24):0.14];
% kD varying:
sIkD = [(10/25):(10/25):10];

%% Kinetic parameters:
k.gD = 0;   % Dilution
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
k.bM = (5/12)*0.30;
k.gM = 1.5;
k.KM = 1;
k.KF = 1; 

% Perturbation
k.P  = [0;      % Time points
        1.05];   % Multiplicative perturbation
k.Pn = 'bC';    % Variable to be perturbed

% Delay
k.do_PID_delay   = 1;
k.tau_PID_delay  = 40;	% minutes
k.N_PID_delay    = 4;   % steps of tau_PID_delay/minutes
k.beta_PID_delay = k.N_PID_delay/k.tau_PID_delay; % rate per step
k.do_PID_delay_feedback = 1;
k.feedback_gain_process = 0.2*k.do_PID_delay_feedback;

options = odeset('AbsTol',1e-6,'RelTol',1e-6);

%% Sweaping:
AT = zeros(length(sIkP),length(sIkD),2) + Inf;
mE = zeros(length(sIkP),length(sIkD),2) - Inf;
mY = zeros(length(sIkP),length(sIkD),6+k.N_PID_delay) - Inf;
for i = 1:length(sIkP)
    for j = 1:length(sIkD)
        [i,j]
        k.bP = sIkP(i)*(4*k.mu);
        k.bD = sIkD(j)*power((1/k.th)*(k.gA/(k.bA*k.gM)),-1);

        % Aproximated steady states
        XC = k.mu*k.Y/k.th;
        X1 = XC*k.gC/k.bC;
        A  = k.bM*k.Y/k.gM;
        M  = ((k.gA*XC)+(k.gA0*A))/k.bA;
        Z1 = ((k.g1*X1)-(k.bP*k.Y*k.mu*k.Y/((k.th*XC)...
            +(k.mu*k.Y)))-(k.bD*A))/k.bI;
        Z2 = k.th*XC/(k.et*Z1);
        % Numerical steady states
        O = k;
        k.P  = [0;1];
        save('Par_ODE.mat','k');
        if(k.do_PID_delay == 0)  
            y_in = [Z1 Z2 X1 XC A M];
        elseif((k.do_PID_delay == 1)&&(k.N_PID_delay == 1)) 
            y_in = [Z1 Z2 X1 XC A M XC];
        elseif((k.do_PID_delay == 1)&&(k.N_PID_delay == 2)) 
            y_in = [Z1 Z2 X1 XC A M XC XC];
        elseif((k.do_PID_delay == 1)&&(k.N_PID_delay == 3))
            y_in = [Z1 Z2 X1 XC A M XC XC XC];
        elseif((k.do_PID_delay == 1)&&(k.N_PID_delay == 4))
            y_in = [Z1 Z2 X1 XC A M XC XC XC XC];
        end
        y_in([y_in<0]) = 1e-3;
        y_in([y_in>10000]) = 10000;
        [t,ySS] = ode15s(@FN_ODE_PID,[0:T],y_in,options);
        k = O;
        SS = ySS(length(t),:);
        if(sum(SS<0)>0)
            'ERROR : Negative numbers in SS'
        end
        clear t ySS Z1 Z2 X1 XC A M y_in

        % ODE - Full system
        save('Par_ODE.mat','k');
        [tF,yF] = ode15s(@FN_ODE_PID,[0:T],SS,options);
        if(sum(sum(yF<0))>0)
            'ERROR : Negative numbers in yF'
        end
        XC  = yF(:,4);
        mY(i,j,:) = mean(yF);
        mE(i,j,:) = [min(XC-mY(i,j,4)) max(XC-mY(i,j,4))];
        
        % Adaptation performance:
        eX = abs(XC-mY(i,j,4))/mY(i,j,4);
        AT(i,j,1) = sum([cumsum([eX([T:-1:0]+1)<=eps]) == [1:(T+1)]']);
        eX = abs(XC-mY(i,j,4))/max(abs(mE(i,j,:)));
        AT(i,j,2) = sum([cumsum([eX([T:-1:0]+1)<=eps]) == [1:(T+1)]']);
        clear tF yF SS XC eX
    end
    save(cat(2,'TEMP_AT_PIDn_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))
end
clear i j O ans
AT = T + 1 - AT;
save(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))
delete(cat(2,'TEMP_AT_PIDn_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))

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
    
    imagesc(sIkD*k.th*k.bA*k.gM/k.gA,sIkP*4*k.mu,at,[0 1000])
        xlabel('\beta_D','FontSize',12)
        ylabel('\beta_P','FontSize',12)
        set(gca,'YDir','normal','FontSize',12)
        colorbar;
clear a b c

%% END