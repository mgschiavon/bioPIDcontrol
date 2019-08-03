%%       FIGURE : PID- control        %%
%%%%%%%% Stochastic simulation %%%%%%%%%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
% Mariana Gómez-Schiavon
% November, 2017

clear;
vol = 2.9e-15;
avo = 6.022e23;
nMM = 1e-9;

% Kinetic parameters:
k.gD = 0;                   % Dilution [1/min]
k.mu = 1;                   % [1/min]
k.Y  = 300*(vol*avo*nMM);   % [molecules]
k.et = 0.01/(vol*avo*nMM);	% [1/(min*molecules)]
k.th = 0.3;                 % [1/min]
k.g1 = 0.1;                 % [1/min]
k.bC = 0.1;                 % [1/min]
k.gC = 0.1;                 % [1/min]
k.bA = 1.5;                 % [1/min]
k.gA = 1.5;                 % [1/min]
k.KA = 1*(vol*avo*nMM);     % [molecules]
k.gA0= 0.1;                 % [1/min]
k.bM = 5/12;                % [1/min]
k.gM = 1.5;                 % [1/min]
k.KM = 1*(vol*avo*nMM);     % [molecules]

% DEFAULT - I-Control
k.bI = 0.06;                % [1/min]
k.bP = 0;                   % [1/min]
k.bD = 0;                   % [1/min]

% Perturbation
k.Pn = 'Y';                 % Variable to be perturbed
k.P  = [0,1000,3000;        % Time points
        1,2,0.2];           % Multiplicative perturbation
% k.Pn = 'bC';                % Variable to be perturbed
% k.P  = [0,1000,3000;        % Time points
%         1,1.5,2];           % Multiplicative perturbation

%% FIGURE
% k.bP = 0.3;       % Proportional control
k.bD = 0.56;       % Derivative control

T = zeros(max(k.P(1,:))+2000,1);
X = zeros(max(k.P(1,:))+2000,6);

t = 0.001;
x = zeros(1,6);

while(t<(max(k.P(1,:))+2000))
    T(ceil(t))   = t;
    X(ceil(t),:) = x;
    [t,x] = FN_SSA_PID(t,x,k);
end
clear t x
save SS_PID_IDc.mat

% Figure:
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 4 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.2 0.15 0.725 0.825])
hold on;
plot(T,X(:,4),'LineWidth',2)
    ylabel('X_C [molecules]','FontSize',14)
    xlabel('Minutes','FontSize',14)
    legend('show')
    ylim([0 3100*(vol*avo*nMM)])
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
else
    plot([T(1),T(length(T))],[0 0]+(k.mu*k.Y/k.th),...
        'LineStyle',':','LineWidth',2,'Color',[0 0 0]+0.7)
end
        xlim([750 5000])
        set(gca,'XTick',[1000:2000:5000],'XTickLabel',[0:2000:5000])
        
%% END