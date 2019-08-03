%%       FIGURE : PID- control        %%
%%%%%%%%%  Adaptation dynamics %%%%%%%%%
%%%%%%% Quantifying performance %%%%%%%%
%%%%%%%%%%  Adaptation Error %%%%%%%%%%%
% Mariana Gómez-Schiavon
% December, 2017

clear;

%% Quantification parameters:
k.Y  = 300;
k.bI = 0.06;
load(cat(2,'AT_PID_Sweap_Y',num2str(k.Y),'_bI',num2str(k.bI),'.mat'))

%% Sweaping:
X2i = k.mu*k.Y/k.th;
Ae = zeros(sweaI1,sweaI2);
for i = 1:sweaI1
    for j = 1:sweaI2
        Ae(i,j) = mY(i,j,4) - X2i;
    end
end

%% Heat-map figure
fig = figure();
fig.Units = 'inches';
fig.Position = [5 2 3 2];
fig.PaperPosition = fig.Position;
C = colormap('jet');
C(1,:)  = [1 1 1];
C(64,:) = [0 0 0];
fig.Colormap = C;    
    imagesc([1:sweaI2]*sweaS2*k.gA/(k.bA*k.gM),...
        [1:sweaI1]*sweaS1*k.th/(4*k.mu),Ae,[-1 1]*(0.1*X2i))
        xlabel('k_D \theta','FontSize',12)
        ylabel('k_P \theta','FontSize',12)
        set(gca,'YDir','normal','FontSize',12,...
            'XTick',[0.2:0.2:0.6])
        colorbar;
clear a b c

%% END