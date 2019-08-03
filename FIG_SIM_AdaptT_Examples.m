%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%% Quantifying performance %%%%%%%%
%%%%%%%%%%%%%%  Examples  %%%%%%%%%%%%%%
% Mariana Gómez-Schiavon
% December, 2017

clear;

%% Quantification parameters:
k.Y  = 300;
k.bI = 0.06;
load(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))

%%
C = colormap('jet');
C(64,:) = [0 0 0];

O = k;
for i = [5,15]
    for j = [5,15]
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
        [t,y] = ode15s(@FN_ODE_PID,[0:T],SS);
        
        at = AT(i,j,2);
        c  = round(64*at/300);
        if(c>64)
            c = 64;
        end
        
        fig = figure();
        fig.Units = 'inches';
        fig.Position = [5 2 2 1.5];
        fig.PaperPosition = fig.Position;
        axes('Position',[0.01 0.01 0.98 0.98])
        hold on;
            plot(t,y(:,4),'LineWidth',2,'Color',C(c,:))
            plot([t(1),t(length(t))],[0 0]+(k.mu*k.Y/k.th),...
                'LineStyle',':','LineWidth',2,'Color',[0 0 0]+0.7)
                xlim([0 T])
                ylim([0.7 1.3]*(k.mu*k.Y/k.th))
                set(gca,'XTick',[0:200:T],...
                    'YTick',[0.8:0.2:1.2]*(k.mu*k.Y/k.th))
                box('on')
                grid('on')
        
        clear t y SS
    end
end
k = O;
clear i j O ans

%%
O = k;
for i = [1,13,sweaI1]
    for j = [1,13,sweaI2]
        k.(sweaP1) = i*sweaS1;
        k.(sweaP2) = j*sweaS2;
        [j*sweaS2*k.gA/(k.bA*k.gM),i*sweaS1*k.th/(4*k.mu)]
    end
end
k = O;

%% END