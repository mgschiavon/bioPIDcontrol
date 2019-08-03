%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%%%%    Showing A~dXC/dt  %%%%%%%%%
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
k.bC = 0.1;
k.gC = 0.1;
k.bA = 1.5;
k.gA = 1.5;
k.KA = 1;
k.gA0= 0.1;
k.bM = 5/12;
k.gM = 1.5;
k.KM = 1;

% DEFAULT - I-Control
k.bI = 0.06;
k.bP = 0;
k.bD = 0;

% Perturbation
% k.Pn = 'Y';             % Variable to be perturbed
% k.P  = [0,1000,3000;    % Time points
%         1,2,0.2];       % Multiplicative perturbation
k.Pn = 'bC';            % Variable to be perturbed
k.P  = [0,1000,3000;    % Time points
        1,1.5,2];       % Multiplicative perturbation

%% SIMULATION
% k.bP = 0.3;       % Proportional control
k.bD = 0.6;       % Derivative control

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
[t,y] = ode45(@FN_ODE_PID,[0 max(k.P(1,:))+2000],[Z1 Z2 X1 XC A M]);

%% Figure
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 8 4];
fig.PaperPosition = fig.Position;
    left_color = [0 0 0];
    right_color = [0.5 0.5 0.5];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    axes(fig,'Position',[0.07 0.15 0.855 0.825])
    
% axes(fig,'Position',[0.1 0.15 0.85 0.85])
hold on;
yyaxis left
    plot(t,-[0;diff(y(:,4))./diff(t)],'LineWidth',2)
        ylabel('~ -dX_C/dt','FontSize',14)
%         xlim([750 4000])
%         ylim([-15.5 4.2])
        set(gca,'YTick','','FontSize',14)
        box('on')
yyaxis right
    plot(t,(y(:,5)-mean(y(:,5)))/mean(y(:,5)),...
        'LineWidth',2,'LineStyle',':')
        ylabel('A','FontSize',14)
        xlabel('Minutes','FontSize',14)
%         xlim([750 4000])
%         ylim([-0.042 0.135])
        set(gca,'YTick','','FontSize',14,'XGrid','on')
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
    
%% Additional format
if(strcmp(k.Pn,'Y'))
    plot([k.P(1,1) k.P(1,2) k.P(1,2) k.P(1,3) k.P(1,3) k.P(1,3)+2000],...
        (k.mu*k.Y/k.th)*[k.P(2,1) k.P(2,1) k.P(2,2) k.P(2,2) k.P(2,3) k.P(2,3)],...
        'LineWidth',2,'LineStyle',':','Color',[0 0 0]+0.7)
end
        xlim([750 5000])
        set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])
        
%% END