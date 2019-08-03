%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%% Quantifying performance %%%%%%%%
%%%%%%%%%%%%%%  Examples  %%%%%%%%%%%%%%
%%%%%%%%% Complex process design %%%%%%%%%
% Mariana Gómez-Schiavon
% June, 2018

clear;

%% Quantification parameters:
k.Y  = 300;
k.bI = 0.02;
load(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))

%%
C = colormap('jet');
C(64,:) = [0 0 0];

O = k;
for i = 10%[1,8,16,24]
    for j = [1,24]%[10]
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
        [t,y] = ode15s(@FN_ODE_PID,[0:T],SS,options);
        if(sum(sum(y<0))>0)
            'ERROR : Negative numbers in yF'
        end
        
        at = AT(i,j,2);
        c  = round(64*at/1000);
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
                ylim([0.95 1.05]*(k.mu*k.Y/k.th))
                set(gca,'XTick',[0:200:T],...
                    'YTick',[0.8:0.05:1.2]*(k.mu*k.Y/k.th))
                xlim([0 T/5])
                box('on')
                grid('on')
        
        clear t y SS
    end
end
k = O;
clear i j O ans

%%
O = k;
for i = 10%[1,8,16,24]
    for j = [1,24]%[10]
        [sIkP(i)*(4*k.mu),sIkD(j)*power((1/k.th)*(k.gA/(k.bA*k.gM)),-1)]
    end
end
k = O;

%% END