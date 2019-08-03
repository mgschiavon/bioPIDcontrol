%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%% Linearized *textbook* system %%%%%
% Mariana Gómez-Schiavon
% December, 2017

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
k.bP = 0.15;
k.bD = 0.3;

% Perturbation
% k.Pn = 'Y';             % Variable to be perturbed
% k.P  = [0,3000,5000;    % Time points
%         1,1.2,0.8];       % Multiplicative perturbation (or [1,2,0.2])
k.Pn = 'bC';            % Variable to be perturbed
k.P  = [0,3000,5000;    % Time points
        1,1.1,1.2];       % Multiplicative perturbation (or [1,1.5,2])

    
%% Find steady states:
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
[t,ySS] = ode45(@FN_ODE_PID,[0:5000],[Z1 Z2 X1 XC A M]);
k = O;
k.SS = ySS(length(t),:);
clear t O ySS

%% ODE - Linearized Textbook
save('Par_ODE.mat','k');
[tX,yX] = ode15s(@FN_ODE_PID_LinTX,[0 max(k.P(1,:))+2000],[Z1 Z2 X1 XC A M 0]);

%% ODE - Linearized
save('Par_ODE.mat','k');
[tL,yL] = ode15s(@FN_ODE_PID_Lin,[0 max(k.P(1,:))+2000],[Z1 Z2 X1 XC A M]);

%% ODE - Full system
save('Par_ODE.mat','k');
[tF,yF] = ode15s(@FN_ODE_PID,[0 max(k.P(1,:))+2000],[Z1 Z2 X1 XC A M]);

%% FIGURE
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 8 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.125 0.15 0.825 0.825])
hold on;
plot(tF,yF(:,4),'LineWidth',2)
plot(tL,k.SS(4)+yL(:,4),'LineWidth',2)
plot(tX,k.SS(4)+yX(:,4),'LineWidth',2)
if(strcmp(k.Pn,'Y'))
    plot([k.P(1,1) k.P(1,2) k.P(1,2) k.P(1,3) k.P(1,3) k.P(1,3)+2000],...
        (k.mu*k.Y/k.th)*[k.P(2,1) k.P(2,1) k.P(2,2) k.P(2,2) k.P(2,3) k.P(2,3)],...
        'LineWidth',2,'LineStyle',':','Color',[0 0 0]+0.7)
else
	plot([tF(1),tF(length(tF))],[0 0]+(k.mu*k.Y/k.th),...
		'LineStyle',':','LineWidth',2,'Color',[0 0 0]+0.7)
end
    legend({'Full system','Linearized','Linearized with Eq.(12)'})
    ylabel('X_C [nM]','FontSize',14)
    xlabel('Minutes','FontSize',14)
    legend('show')
	set(gca,'XTick',[1000:2000:5000]+2000,'XTickLabel',[0:2000:5000])
	xlim([750 5000]+2000)
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
    annotation('textbox',[0.18 0.7 0.45 0.2],...
        'String',{cat(2,'\mu=',num2str(k.mu),...
        '; Y=',num2str(k.Y),...
        '; \eta=',num2str(k.et),...
        '; \gamma_D=',num2str(k.gD),...
        '; \theta=',num2str(k.th),...
        '; \gamma_1=',num2str(k.g1),...
        '; \beta_C=',num2str(k.bC),...
        '; \gamma_C=',num2str(k.gC),...
        '; \beta_A=',num2str(k.bA),...
        '; \gamma_A=',num2str(k.gA),...
        '; K_A=',num2str(k.KA),...
        '; \gamma_{A_0}=',num2str(k.gA0),...
        '; \beta_MY=',num2str(k.bM*k.Y),...
        '; \gamma_M=',num2str(k.gM),...
        '; K_M=',num2str(k.KM),...
        '; \beta_I=',num2str(k.bI),...
        '; \beta_P=',num2str(k.bP),...
        '; \beta_D=',num2str(k.bD))},...
        'FitBoxToText','off');

%% END