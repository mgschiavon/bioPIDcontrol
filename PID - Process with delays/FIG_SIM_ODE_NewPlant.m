%%       FIGURE : PID-control         %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%%%% Complex process design %%%%%%%%%
% Mariana Gómez-Schiavon
% June, 2018

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
k.bM = (5/12)*0.30;
k.gM = 1.5;
k.KM = 1;
k.KF = 1; 

% DEFAULT - I-Control
k.bI = 0.06;
k.bP = 0.0;     % I
k.bD = 0.0;     % I 
%k.bP = 0.15;	% PID
%k.bD = 0.28;	% PID
%k.bP = 0.3;	% PI
%k.bD = 0.0;	% PI
%k.bP = 0.0;	% ID
%k.bD = 0.56;	% ID

k.Y_shift_D_linear = 1; 	% 1 - keep normal, 0 - no Y shift in the D control

k.do_PID_delay  = 1;
k.tau_PID_delay = 40;	% minutes
k.N_PID_delay   = 4;    % steps of tau_PID_delay/minutes
if (k.do_PID_delay == 0)
	k.N_PID_delay   = 0;
    k.tau_PID_delay = 0;
end
k.beta_PID_delay = k.N_PID_delay/k.tau_PID_delay; % rate per step

k.do_PID_delay_feedback = 1;
k.feedback_gain_process = 0.2*k.do_PID_delay_feedback;

k.kP_test = 0.03;
k.kI_test = 0.02;
k.kD_test = 4;

k.bI   = k.kI_test;
k.bP   = (4*k.mu)*k.kP_test;
k.bD   = k.kD_test*power((1/k.th)*(k.gA/(k.bA*k.gM)),-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Mariana:  Start sweeping about his point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bI_paper = k.bI
bP_paper = k.bP
bD_paper = k.bD


do_ode15s = 1;      % 1 - Yes, 0 - NO
t_extra   = 0*2000; % extra ODE time
k.Pn = 'bC';            % Variable to be perturbed
k.P  = [0,3000;    % Time points
        1,1.05];   % Multiplicative perturbation (or [1,1.5,2])
    
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
if (k.do_PID_delay == 0)  
    y_in = [Z1 Z2 X1 XC A M];
elseif (k.do_PID_delay == 1)&(k.N_PID_delay == 1)  
    y_in = [Z1 Z2 X1 XC A M XC];
elseif (k.do_PID_delay == 1)&(k.N_PID_delay == 2)  
    y_in = [Z1 Z2 X1 XC A M XC XC];
elseif (k.do_PID_delay == 1)&(k.N_PID_delay == 3)
    y_in = [Z1 Z2 X1 XC A M XC XC XC];
elseif (k.do_PID_delay == 1)&(k.N_PID_delay == 4)
    y_in = [Z1 Z2 X1 XC A M XC XC XC XC];
end;
clear Z1 Z2 X1 XC A M
options = odeset('AbsTol',1e-5,'RelTol',1e-5);
if (do_ode15s == 0)    
    [t,ySS] = ode45(@FN_ODE_PID,[0:5000],y_in);
elseif (do_ode15s == 1)        
    [t,ySS] = ode15s(@FN_ODE_PID,[0:5000],y_in,options);
end;
k = O;
k.SS = ySS(length(t),:);
clear t O ySS

%     Z1s = k.SS(1);
%     Z2s = k.SS(2);
%     kI  = k.bI*Z1s/(Z1s+Z2s)
%     kP  = k.bP/(4*k.mu)
%     kD  = (1/k.th)*(k.bD*k.gA/(k.bA*k.gM))

%% ODE - Full system
save('Par_ODE.mat','k');
if (do_ode15s == 0)    
    [tF,yF] = ode45(@FN_ODE_PID,[0 max(k.P(1,:))+2000+t_extra],k.SS);
elseif (do_ode15s == 1)        
    [tF,yF] = ode15s(@FN_ODE_PID,[0 max(k.P(1,:))+2000+t_extra],k.SS,options);
end

%% FIGURE
index_label = {'Z_1','Z_2','X_1','X_C','A','M'};
index_plot  = 4;	%[1-Z1,2-Z2,3-X1,4-XC,5-A,6-M]
fig = figure();
fig.Units = 'inches';
fig.Position = [8 2 8 4];
fig.PaperPosition = fig.Position;
axes(fig,'Position',[0.125 0.175 0.825 0.75])
hold on;
    if (index_plot==4)
        if(strcmp(k.Pn,'Y'))
            plot([k.P(1,1) k.P(1,2) k.P(1,2) k.P(1,3) k.P(1,3) k.P(1,3)+2000],...
                (k.mu*k.Y/k.th)*[k.P(2,1) k.P(2,1) k.P(2,2) k.P(2,2) k.P(2,3) k.P(2,3)],...
                'LineWidth',2,'LineStyle',':','Color',[0 0 0]+0.7)
        else 
            plot([tF(1),tF(length(tF))],[0 0]+(k.mu*k.Y/k.th),... 
                'LineStyle',':','LineWidth',2,'Color',[0 0 0]+0.7)
        end
        legend('Set-point: \muY/\theta')
    end
    plot(tF,yF(:,index_plot),'LineWidth',2)
hold off;
xlabel('Minutes')
ylabel(cat(2,index_label{index_plot},' [nM]'),'FontSize',14)
set(gca,'XTick',[1000:1000:(max(k.P(1,:))+t_extra)]+2000,...
    'XTickLabel',[0:1000:(max(k.P(1,:))+t_extra)])
    xlim([750 (max(k.P(1,:))+t_extra)]+2000)
    if (k.do_PID_delay == 0)
        title(strcat('Process delay: 0 minutes'));
    elseif (k.do_PID_delay == 1)
        title(cat(2,'Process delay: ',num2str(k.tau_PID_delay),...
            ' minutes in ',num2str(k.N_PID_delay),' steps'));
    end;
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
        '; \beta_D=',num2str(k.bD),...
        '; k_I=',num2str(k.kI_test),...
        '; k_P=',num2str(k.kP_test),...
        '; k_D=',num2str(k.kD_test))},...
        'FitBoxToText','off');

%% END
