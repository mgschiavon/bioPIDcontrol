%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%%  All species & Z2 inset %%%%%%%%
%%%%%%% Effect of alpha_P not 1 %%%%%%%%
% Mariana Gómez-Schiavon
% November, 2017

clear;
% Kinetic parameters:
k.gD = 0;   % Dilution
k.mu = 1;
k.Y  = 300;
k.et = 0.01;
k.th = 0.3;
k.g1 = 0.1;
k.bC = 0.2;
k.gC = 0.1;
k.bA = 1.5;
k.gA = 1.5;
k.KA = 1;
k.gA0= 0.1;
k.bM = 5/12;
k.gM = 1.5;
k.KM = 1;

k.aP  = 1.5;

% DEFAULT - I-Control
k.bI = 0.06;
k.bP = 0;
k.bD = 0;

% Perturbation -- NONE
k.Pn = 'Y';       % Variable to be perturbed
k.P  = [0;        % Time points
        1];       % Multiplicative perturbation

%% Simulation
k.bP = 0.3;       % Proportional control
% k.bD = 1.5;       % Derivative control

% Aproximated steady states
XC = k.mu*k.Y/k.th;
X1 = XC*k.gC/k.bC;
A  = k.bM*k.Y/k.gM;
M  = ((k.gA*XC)+(k.gA0*A))/k.bA;
Z1 = ((k.g1*X1)-(k.bP*k.Y*k.mu*k.Y/((k.th*XC)...
    +(k.mu*k.Y)))-(k.bD*A))/k.bI;
Z2 = k.th*XC/(k.et*Z1);

% ODE
save('Par_ODE.mat','k');
[t,y] = ode45(@FN_ODE_PID,[0 400],abs([Z1 Z2 X1 XC A M])+rand(1,6));

%% Figure:
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.825])
hold on;

plot(t,y,'LineWidth',2)    
    ylabel('Concentration [nM]','FontSize',14)
    xlabel('Minutes','FontSize',14)
%     legend({'Z_1','Z_2','X_1','X_C','A','M'})
    ylim([0 3100])
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
%     annotation('textbox',[0.18 0.7 0.45 0.2],...
%         'String',{cat(2,'\mu=',num2str(k.mu),...
%         '; Y=',num2str(k.Y),...
%         '; \eta=',num2str(k.et),...
%         '; \gamma_D=',num2str(k.gD),...
%         '; \theta=',num2str(k.th),...
%         '; \gamma_1=',num2str(k.g1),...
%         '; \beta_C=',num2str(k.bC),...
%         '; \gamma_C=',num2str(k.gC),...
%         '; \beta_A=',num2str(k.bA),...
%         '; \gamma_A=',num2str(k.gA),...
%         '; K_A=',num2str(k.KA),...
%         '; \gamma_{A_0}=',num2str(k.gA0),...
%         '; \beta_MY=',num2str(k.bM*k.Y),...
%         '; \gamma_M=',num2str(k.gM),...
%         '; K_M=',num2str(k.KM),...
%         '; \beta_I=',num2str(k.bI),...
%         '; \beta_P=',num2str(k.bP),...
%         '; \beta_D=',num2str(k.bD))},...
%         'FitBoxToText','off');

axes(fig,'Position',[0.6 0.65 0.3 0.3])
    plot(t,y(:,2),'LineWidth',2)    
        ylabel('Z_2 [nM]','FontSize',14)
        xlabel('Minutes','FontSize',14)
        box('on')
        set(gca,'FontSize',14,'XGrid','on')

%% END