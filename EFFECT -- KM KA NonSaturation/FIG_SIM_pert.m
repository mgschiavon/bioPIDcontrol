%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
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

%% FIGURE
% k.bP = 0.3;       % Proportional control
k.bD = 0.56;       % Derivative control

% Figure:
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.825])
hold on;

C = colormap('cool');
pI = [100];
for i = 1:length(pI)
    k.KM  = pI(i);
    
    % Aproximated steady states
    XC = k.mu*k.Y/k.th;
    X1 = XC*k.gC/k.bC;
    A  = k.bM*k.Y/k.gM;
    M  = ((k.gA*XC)+(k.gA0*A))/k.bA;
    Z1 = ((k.g1*X1)-(k.bP*k.Y*k.mu*k.Y/((k.th*XC)...
        +(k.mu*k.Y)))-(k.bD*A))/k.bI;
    Z2 = k.th*XC/(k.et*Z1);
    
    % Correcting effect of large KM by adjusting gM:
    % NOTE: Assumming (k.gA0/k.bA) is small enough...
%     k.gM = k.gM/((k.gA*XC/k.bA)/((k.gA*XC/k.bA)+k.KM));

    % ODE
    save('Par_ODE.mat','k');
    [t,y] = ode45(@FN_ODE_PID,[0 max(k.P(1,:))+2000],abs([Z1 Z2 X1 XC A M]));
    
    Y(i).t = t;
    Y(i).y = y;

    plot(Y(i).t,Y(i).y(:,4),'LineWidth',1,...
        'Color',C(i*floor(64/length(pI)),:),...
        'DisplayName',num2str(k.KM))
end
    ylabel('X_C [nM]','FontSize',14)
    xlabel('Minutes','FontSize',14)
    legend('show')
    ylim([0 3100])
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
    
%% Additional format
if(strcmp(k.Pn,'Y'))
    plot([k.P(1,1) k.P(1,2) k.P(1,2) k.P(1,3) k.P(1,3) k.P(1,3)+2000],...
        (k.mu*k.Y/k.th)*[k.P(2,1) k.P(2,1) k.P(2,2) k.P(2,2) k.P(2,3) k.P(2,3)],...
        'LineWidth',2,'LineStyle',':','Color',[0 0 0]+0.7)
end
        xlim([750 5000])
        set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])

%%
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.825])
hold on;
for i = 1:length(pI)
    k.KM  = pI(i);
    plot(Y(i).t,Y(i).y,'LineWidth',2)
end
    ylabel('[nM]','FontSize',14)
    xlabel('Minutes','FontSize',14)
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
    xlim([750 5000])
    set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])
    legend({'Z_1','Z_2','X_1','X_C','A','M'})

%%
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.225])
hold on;

for i = 1:length(pI)
    k.KM  = pI(i);
    hold on;
    plot(Y(i).t,k.KM./Y(i).y(:,6),'LineWidth',1,...
        'Color',C(i*floor(64/length(pI)),:),...
        'DisplayName',num2str(k.KM))
end
    ylabel('K_M/M ','FontSize',14)
    xlabel('Minutes','FontSize',14)
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
    xlim([750 5000])
    set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])

axes(fig,'Position',[0.2 0.65 0.725 0.225])
hold on;

for i = 1:length(pI)
    k.KM  = pI(i);
    hold on;
    plot(Y(i).t,Y(i).y(:,6),'LineWidth',1,...
        'Color',C(i*floor(64/length(pI)),:),...
        'DisplayName',num2str(k.KM))
end
    ylabel('M ','FontSize',14)
    xlabel('Minutes','FontSize',14)
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
    xlim([750 5000])
    set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])

%% END