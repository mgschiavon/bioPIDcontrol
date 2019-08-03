%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%% Simulate a perturbation %%%%%%%%
%%%%%%%% New controller design %%%%%%%%%
% Mariana Gómez-Schiavon
% December, 2017

fac_bZ = 2;
fac_bI = 0.96;

% Kinetic parameters:
k.g   = 0;      % Dilution
k.Y   = 300;
k.g1  = 0.1;
k.mu  = 1;
k.th  = 0.3;
k.bA  = 1.5;
k.gA  = 1.5;
k.KA  = 1;
k.gA0 = 0.1;
k.bM  = 0.4167;
k.gM  = 1.5;
k.KM  = 1;
k.bZ  = fac_bZ*2;
k.gZ  = 1;
k.aZ  = fac_bZ*1;
k.KZ  = 1.0;
k.bC  = 0.1;
k.gC  = 0.1;


% DEFAULT - I-Control
k.bI  = fac_bI*0.12/fac_bZ;
k.bP = 0;
k.bD = 0;

% Perturbation
k.Pn = 'Y';             % Variable to be perturbed
k.P  = [0,3000,5000;    % Time points
        1,2,0.2];       % Multiplicative perturbation
% k.Pn = 'bC';            % Variable to be perturbed
% k.P  = [0,3000,5000;    % Time points
%         1,1.5,2];       % Multiplicative perturbation

%% FIGURE
k.bP = 0.15;       % Proportional control
k.bD = 0.3;       % Derivative control

% Figure:
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.825])
hold on;

C = colormap('cool');
if(strcmp(k.Pn,'Y'))
    pI = k.Y;
else
    pI = [60,180,300,540];%[60:120:600];
end
for i = 1:length(pI)
    k.Y  = pI(i);
    
    % ODE
    save('Par_ODE.mat','k');
    [t,y] = ode45(@FN_ODE_PID_TFC,[0 max(k.P(1,:))+2000],zeros(1,5));
    
    Y(i).t = t;
    Y(i).y = y;

    plot(Y(i).t,Y(i).y(:,5),'LineWidth',1,...
        'Color',C(i*floor(64/length(pI)),:),...
        'DisplayName',num2str(k.mu*k.Y))
end
    ylabel('X_C [nM]','FontSize',14)
    xlabel('Minutes','FontSize',14)
    legend('show')
    ylim([0 3100])
    box('on')
    set(gca,'FontSize',14,'XGrid','on')
    annotation('textbox',[0.18 0.7 0.45 0.2],...
        'String',{cat(2,'\gamma=',num2str(k.g),...
        '; Y=',num2str(k.Y),...
        '; \gamma_1=',num2str(k.g1),...
        '; \mu=',num2str(k.mu),...
        '; \theta=',num2str(k.th),...
        '; \beta_A=',num2str(k.bA),...
        '; \gamma_A=',num2str(k.gA),...
        '; K_A=',num2str(k.KA),...
        '; \gamma_{A_0}=',num2str(k.gA0),...
        '; \beta_M=',num2str(k.bM),...
        '; \gamma_M=',num2str(k.gM),...
        '; K_M=',num2str(k.KM),...
        '; \beta_Z=',num2str(k.bZ),...
        '; \gamma_Z=',num2str(k.gZ),...
        '; \alpha_Z=',num2str(k.aZ),...
        '; K_Z=',num2str(k.KZ),...
        '; \beta_C=',num2str(k.bC),...
        '; \gamma_C=',num2str(k.gC),...
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
        xlim([750 5000]+2000)
        set(gca,'XTick',[1000:2000:5000]+2000,'XTickLabel',[0:2000:5000])
        
%% END